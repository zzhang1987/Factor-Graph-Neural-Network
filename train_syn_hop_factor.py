import lib
import torch
import argparse
import sys
import utils
import logging
import datetime
from lib.model.mpnn import factor_mpnn
import numpy as np
import time
import os
from tensorboardX import SummaryWriter
from utils.types import str2bool
from tqdm import tqdm
import statistics as st


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--chain_length',
                        type=int,
                        default=30,
                        help="The length of generated chain structured MRF")
    parser.add_argument(
        '--hop_cap',
        type=int,
        default=5,
        help="The seed to generate parameter of budget factors")
    parser.add_argument('--nfactors',
                        type=int,
                        default=8,
                        help="Number of higher order factors")
    parser.add_argument('--hop_order',
                        type=int,
                        default=9,
                        help="Order of higher order factors")
    parser.add_argument('--train_epoches',
                        type=int,
                        default=10,
                        help="training epoches")
    parser.add_argument('--model_path',
                        type=str,
                        default='gnn',
                        help="Saved model path")
    parser.add_argument('--model_name',
                        type=str,
                        default='mp_nn_factor',
                        help="model name (PointNet, GNN)")
    parser.add_argument('--neighbour',
                        type=int,
                        default=9,
                        help="number of neighbour in the graph")

    parser.add_argument('--log_level',
                        type=str,
                        default='info',
                        help="log level")

    parser.add_argument('--use_cuda',
                        type=str2bool,
                        default=True,
                        help="Use cuda or not")

    parser.add_argument('--train_path',
                        type=str,
                        default="synthetic_data/hop_train.dat",
                        help="path of the training dataset")

    parser.add_argument('--test_path',
                        type=str,
                        default="synthetic_data/hop_test.dat",
                        help="path of the testing dataset")

    parser.add_argument('--train_size',
                        type=int,
                        default=90000,
                        help="size of training dataset")

    parser.add_argument('--test_size',
                        type=int,
                        default=10000,
                        help="size of testing dataset")

    parser.add_argument('--batch_size', type=int, default=32)

    return parser.parse_args()


def generate_knn_table(n, k, knn=False):
    if k % 2 == 0:
        k = k + 1
    nn_idx = np.zeros([n, k]).astype(np.int64)
    if knn: efeature = np.zeros([1, n, k]).astype(np.float32)
    hk = k // 2
    for i in range(n):
        for idx, j in enumerate(range(i - hk, i + hk)):
            if j < 0:
                j = 0
            if j >= n:
                j = n - 1
            nn_idx[i, idx] = j
            if knn: efeature[0, i, idx] = i - j
    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    if knn:
        efeature = torch.from_numpy(np.expand_dims(efeature, 0))
        return nn_idx, efeature
    else:
        return nn_idx


def generate_pw_factor_table(n):
    nn_idx = np.zeros([2 * n, 2]).astype(np.int64)
    efeature = np.zeros([3, 2 * n, 2]).astype(np.float32)

    for i in range(n):
        nn = [(i - 1) % n, i]
        for idx, neighbour in enumerate(nn):
            efeature[0, i, idx] = 1
            nn_idx[i, idx] = n + neighbour
            efeature[2, i, idx] = (i - neighbour + 0.5) * 2

        nn = [i, (i + 1) % n]
        for idx, neighbour in enumerate(nn):
            efeature[1, n + i, idx] = 1
            nn_idx[n + i, idx] = neighbour
            efeature[2, n + i, idx] = (i - neighbour + 0.5) * 2

    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    efeature = torch.from_numpy(np.expand_dims(efeature, 0))

    return nn_idx, efeature


def generate_high_factor_table(n, k):
    nn_idx = np.zeros([n << 1, k]).astype(np.int64)
    efeature = np.zeros([2, n << 1, k]).astype(np.float32)

    hk = k >> 1
    for i in range(n):
        for idx in range(k):
            neighbor = (i + idx - hk + n) % n
            efeature[0, i, idx] = 1
            nn_idx[i, idx] = neighbor + n

            efeature[1, n + i, idx] = 1
            nn_idx[n + i, idx] = neighbor

    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    efeature = torch.from_numpy(np.expand_dims(efeature, 0))
    return nn_idx, efeature


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    np.random.seed(((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) +
                   ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24))


def main():
    args = parse_args()
    subdir = f'train_syn_hop_factor_{args.model_name}_nn_{args.neighbour}_at_{datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}'
    utils.init_logger('./logs/', subdir, print_log=False)
    logging.info(str(args))
    logdir = f'./tf_logs/{subdir}'
    writer = SummaryWriter(log_dir=logdir)

    nfeature_dim = 2
    if args.model_name == 'mp_nn_factor':
        model = factor_mpnn(nfeature_dim, [nfeature_dim**2, args.hop_order],
                            [64, 64, 128, 128, 256, 256, 128, 128, 64, 64, 2],
                            [16, 16])

        emodel_pw = torch.nn.Sequential(torch.nn.Conv2d(3, 64, 1),
                                        torch.nn.ReLU(inplace=True),
                                        torch.nn.Conv2d(64, 16, 1))
        emodel_high = torch.nn.Sequential(torch.nn.Conv2d(2, 64, 1),
                                          torch.nn.ReLU(inplace=True),
                                          torch.nn.Conv2d(64, 16, 1))

    def get_model_description():
        return str(model) + str(emodel_pw) + str(emodel_high)

    logging.info('model {} created'.format(get_model_description()))

    cap = args.hop_cap

    nn_idx_pw, efeature_pw = generate_pw_factor_table(args.chain_length)
    nn_idx_high, efeature_high = generate_high_factor_table(
        args.chain_length, args.hop_order)

    if args.use_cuda:
        nn_idx_pw = nn_idx_pw.cuda()
        efeature_pw = efeature_pw.cuda()
        nn_idx_high = nn_idx_high.cuda()
        efeature_high = efeature_high.cuda()

        model.cuda()
        emodel_pw.cuda()
        emodel_high.cuda()

    parameters = list(model.parameters()) + \
                list(emodel_pw.parameters()) + \
                list(emodel_high.parameters())

    # train_data_set = lib.data.RandomPGMHop(args.chain_length,
    #                                       ret_efeature_pw=False)

    # dataloader = torch.utils.data.DataLoader(train_data_set,
    #                                          batch_size=args.batch_size,
    #                                          shuffle=True,
    #                                          num_workers=8,
    #                                          worker_init_fn=worker_init_fn)

    train_dataset = lib.data.RandomPGMData(args.train_path,
                                           pgm_type="hops",
                                           size=args.train_size)
    test_dataset = lib.data.RandomPGMData(args.test_path,
                                          pgm_type="hops",
                                          size=args.test_size)

    train_loader = torch.utils.data.DataLoader(train_dataset,
                                               batch_size=args.batch_size,
                                               shuffle=True,
                                               num_workers=8,
                                               worker_init_fn=worker_init_fn)
    test_loader = torch.utils.data.DataLoader(test_dataset,
                                              batch_size=args.batch_size,
                                              shuffle=True,
                                              num_workers=8,
                                              worker_init_fn=worker_init_fn)

    optimizer = torch.optim.Adam(parameters, lr=3e-3)
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_lambda=lambda x: max(0.98**x, 1e-6))
    start_epoch = 0
    gcnt = 0
    if os.path.exists(args.model_path):
        ckpt = torch.load(args.model_path)
        model.load_state_dict(ckpt['model_state_dict'])
        emodel_pw.load_state_dict(ckpt['emodel_pw_state_dict'])
        emodel_high.load_state_dict(ckpt['emodel_high_state_dict'])
        optimizer.load_state_dict(ckpt['optimizer_state_dict'])
        scheduler.load_state_dict(ckpt['lr_sche'])

        start_epoch = ckpt['epoch']
        gcnt = ckpt['gcnt']

    def get_model_dict():
        return {
            'model_state_dict': model.state_dict(),
            'emodel_pw_state_dict': emodel_pw.state_dict(),
            'emodel_high_state_dict': emodel_high.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'lr_sche': scheduler.state_dict(),
            'epoch': epoch,
            'gcnt': gcnt
        }

    epoch = 0
    for epoch in tqdm(range(start_epoch, args.train_epoches)):
        torch.save(
            get_model_dict(),
            '{}_nn_factor_{}_epoches_{}.pt'.format(args.model_name,
                                                   args.neighbour, epoch))

        logging.info('save train result to {}'.format(
            '{}_nn_factor_{}_epoches_{}.pt'.format(args.model_name,
                                                   args.neighbour, epoch)))
        scheduler.step()

        loss_seq = []
        acc_seq = []
        acc_lp_seq = []
        for bcnt, (nfeature, pws, hops, nlabel,
                   lp_label) in tqdm(enumerate(train_loader)):
            optimizer.zero_grad()
            if args.use_cuda:
                nfeature, pws, hops, nlabel, lp_label \
                    = nfeature.cuda(), pws.cuda(), hops.cuda(), nlabel.cuda(), lp_label.cuda()
            if len(nfeature.shape) == 3:
                nfeature = nfeature.unsqueeze(-1)

            etype_pw = emodel_pw(efeature_pw)
            etype_high = emodel_high(efeature_high)
            bsize = nfeature.shape[0]

            pred, _ = model(nfeature, [pws, hops],
                            [[
                                nn_idx_pw.repeat(bsize, 1, 1),
                                etype_pw.repeat(bsize, 1, 1, 1)
                            ],
                             [
                                 nn_idx_high.repeat(bsize, 1, 1),
                                 etype_high.repeat(bsize, 1, 1, 1)
                             ]])

            pred = pred.squeeze(-1).permute(0, 2, 1).contiguous()
            loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                     nlabel.view(-1))
            loss.backward()
            torch.nn.utils.clip_grad_norm(parameters, 1.0)

            optimizer.step()
            loss_seq.append(loss.item())
            gcnt += 1

            pred_int = pred.argmax(dim=-1)
            all_correct = torch.sum(pred_int == nlabel)
            lp_correct = torch.sum(lp_label == nlabel)
            acc = all_correct.item() / np.prod(nlabel.shape)
            lp_acc = lp_correct.item() / np.prod(nlabel.shape)

            acc_lp_seq.append(lp_acc)
            acc_seq.append(acc)

            if gcnt % 10 == 0:
                logging.info(
                    'epoch = {} bcnt = {} loss = {} acc = {} lp_acc={}'.format(
                        epoch, bcnt, np.mean(loss_seq), np.mean(acc_seq),
                        np.mean(acc_lp_seq)))
                writer.add_scalar('syn_train/loss', loss.item(), gcnt)
                writer.add_scalar('syn_train/acc', acc, gcnt)
                writer.add_scalar('syn_train/lp_acc', lp_acc, gcnt)
                loss_seq = []
                acc_seq = []
                acc_lp_seq = []

    if epoch == args.train_epoches - 1:
        epoch = args.train_epoches
        torch.save(
            get_model_dict(),
            '{}_nn_factor_{}_epoches_{}.pt'.format(args.model_name,
                                                   args.neighbour, epoch))

        logging.info('save train result to {}'.format(
            '{}_nn_factor_{}_epoches_{}.pt'.format(args.model_name,
                                                   args.neighbour, epoch)))
        logging.info('training done!')

    loss_seq = []
    acc_seq = []
    acc_lp_seq = []
    acc_global = []
    acc_lp_global = []
    gcnt = 0
    accum_acc = 0
    accum_acc_lp = 0
    model.eval()
    emodel_high.eval()
    emodel_pw.eval()
    for bcnt, (nfeature, pws, hops, nlabel,
               lp_label) in tqdm(enumerate(test_loader)):
        if args.use_cuda:
            nfeature, pws, hops, nlabel, lp_label \
                = nfeature.cuda(), pws.cuda(), hops.cuda(), nlabel.cuda(), lp_label.cuda()
        if len(nfeature.shape) == 3:
            nfeature = nfeature.unsqueeze(-1)

        etype_pw = emodel_pw(efeature_pw)
        etype_high = emodel_high(efeature_high)
        bsize = nfeature.shape[0]

        pred, _ = model(
            nfeature, [pws, hops],
            [[nn_idx_pw.repeat(bsize, 1, 1),
              etype_pw.repeat(bsize, 1, 1, 1)],
             [
                 nn_idx_high.repeat(bsize, 1, 1),
                 etype_high.repeat(bsize, 1, 1, 1)
             ]])

        pred = pred.squeeze(-1).permute(0, 2, 1).contiguous()
        loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                 nlabel.view(-1))
        torch.nn.utils.clip_grad_norm(parameters, 1.0)

        loss_seq.append(loss.item())
        gcnt += 1

        pred_int = pred.argmax(dim=-1)
        all_correct = torch.sum(pred_int == nlabel)
        lp_correct = torch.sum(lp_label == nlabel)
        acc = all_correct.item() / np.prod(nlabel.shape)
        lp_acc = lp_correct.item() / np.prod(nlabel.shape)
        acc_global.append(acc)
        acc_lp_global.append(lp_acc)

        acc_lp_seq.append(lp_acc)
        acc_seq.append(acc)
        accum_acc += acc
        accum_acc_lp += lp_acc

        if gcnt % 10 == 0:
            logging.info(
                'epoch = {} bcnt = {} loss = {} acc = {} lp_acc={}'.format(
                    epoch, bcnt, np.mean(loss_seq), np.mean(acc_seq),
                    np.mean(acc_lp_seq)))
            writer.add_scalar('syn_test/loss', loss.item(), gcnt)
            writer.add_scalar('syn_test/acc', acc, gcnt)
            writer.add_scalar('syn_test/lp_acc', lp_acc, gcnt)
            loss_seq = []
            acc_seq = []
            acc_lp_seq = []
    logging.info(
        f'testing result: acc = {accum_acc / gcnt}, acc_lp = {accum_acc_lp / gcnt}'
    )
    logging.info(
        f'stddev = {st.stdev(acc_global)}, stddev_lp = {st.stdev(acc_lp_global)}'
    )


if __name__ == '__main__':
    main()
