import lib
import datetime
import torch
import argparse
import sys
import utils
import logging
from tensorboardX import SummaryWriter
from lib.model.mpnn import mp_sequential, mp_conv_residual, mp_conv_type, mp_conv_v2, global_pooling
import numpy as np
import time
import os
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
                        default='mp_nn',
                        help="model name (PointNet, GNN)")
    parser.add_argument('--neighbour',
                        type=int,
                        default=8,
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
                        default="synthetic_data/raw_train.dat",
                        help="path of the training dataset")

    parser.add_argument('--test_path',
                        type=str,
                        default="synthetic_data/raw_test.dat",
                        help="path of the testing dataset")

    parser.add_argument('--train_size', type=int, default=90000,
                        help="size of training dataset")

    parser.add_argument('--test_size', type=int, default=10000,
                        help="size of testing dataset")

    parser.add_argument('--batch_size', type=int, default=32)

    return parser.parse_args()


def generate_knn_table(n, k):
    nn_idx = np.zeros([n, k]).astype(np.int64)
    efeature = np.zeros([1, n, k]).astype(np.float32)
    hk = k // 2
    for i in range(n):
        arr = list(range(i-hk, i)) + list(range(i+1, i+hk))
        for idx, j in enumerate(arr):
            if j < 0:
                j = 0
            if j >= n:
                j = n - 1
            nn_idx[i, idx] = j
            efeature[0, i, idx] = i - j
    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    efeature = torch.from_numpy(np.expand_dims(efeature, 0))
    return nn_idx, efeature


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    np.random.seed(((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) +
                   ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24))


def main():
    args = parse_args()
    subdir = f'raw_nn_{args.neighbour}_at_{datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}'
    utils.init_logger('./logs/', subdir, print_log=False)
    logging.info(str(args))

    writer = SummaryWriter(log_dir=f'./tf_logs/{subdir}')

    nfeature_dim = 2
    print(nfeature_dim)
    if args.model_name == 'mp_nn':
        model = mp_sequential(
            mp_conv_v2(nfeature_dim,
                       64,
                       16,
                       extension=mp_conv_type.ORIG_WITH_NEIGHBOR),
            mp_conv_residual(64, 64, 16), torch.nn.Conv2d(64, 128, 1),
            torch.nn.BatchNorm2d(128), torch.nn.ReLU(inplace=True),
            mp_conv_residual(128, 64, 16), torch.nn.Conv2d(128, 256, 1),
            torch.nn.BatchNorm2d(256), torch.nn.ReLU(inplace=True),
            mp_conv_residual(256, 64, 16), torch.nn.Conv2d(256, 128, 1),
            torch.nn.BatchNorm2d(128), torch.nn.ReLU(inplace=True),
            mp_conv_residual(128, 64, 16), torch.nn.Conv2d(128, 64, 1),
            torch.nn.BatchNorm2d(64), torch.nn.ReLU(inplace=True),
            mp_conv_residual(64, 64, 16), torch.nn.Conv2d(64, 2, 1))
        emodel = torch.nn.Sequential(torch.nn.Conv2d(1, 64, 1),
                                     torch.nn.ReLU(inplace=True),
                                     torch.nn.Conv2d(64, 16, 1))

    elif args.model_name == 'mp_nn_comp':
        model = mp_sequential(
            mp_conv_v2(nfeature_dim,
                       64,
                       16,
                       extension=mp_conv_type.ORIG_WITH_NEIGHBOR),
            mp_conv_residual(64, 64, 16), torch.nn.Conv2d(64, 128, 1),
            torch.nn.BatchNorm2d(128), torch.nn.ReLU(inplace=True),
            mp_conv_residual(128, 64, 16), torch.nn.Conv2d(128, 256, 1),
            torch.nn.BatchNorm2d(256), torch.nn.ReLU(inplace=True),
            mp_conv_residual(256, 64, 16), mp_conv_residual(256, 64, 16),
            mp_conv_residual(256, 64, 16), mp_conv_residual(256, 64, 16),
            mp_conv_residual(256, 64, 16), torch.nn.Conv2d(256, 128, 1),
            torch.nn.BatchNorm2d(128), torch.nn.ReLU(inplace=True),
            mp_conv_residual(128, 64, 16), torch.nn.Conv2d(128, 64, 1),
            torch.nn.BatchNorm2d(64), torch.nn.ReLU(inplace=True),
            mp_conv_residual(64, 64, 16), torch.nn.Conv2d(64, 2, 1))
        emodel = torch.nn.Sequential(torch.nn.Conv2d(1, 64, 1),
                                     torch.nn.ReLU(inplace=True),
                                     torch.nn.Conv2d(64, 16, 1))

    elif args.model_name == 'simple_gnn':
        model = mp_sequential(
            mp_conv_v2(nfeature_dim,
                       64,
                       16,
                       extension=mp_conv_type.ORIG_WITH_NEIGHBOR),
            mp_conv_residual(64, 64, 16), torch.nn.Conv2d(64, 2, 1))
        emodel = torch.nn.Sequential(torch.nn.Conv2d(1, 64, 1),
                                     torch.nn.ReLU(inplace=True),
                                     torch.nn.Conv2d(64, 16, 1))
    elif args.model_name == 'iid':
        model = mp_sequential(torch.nn.Conv2d(nfeature_dim, 64, 1),
                              torch.nn.ReLU(True), torch.nn.Conv2d(64, 2, 1))
        emodel = torch.nn.Sequential(torch.nn.Conv2d(1, 64, 1),
                                     torch.nn.ReLU(inplace=True),
                                     torch.nn.Conv2d(64, 16, 1))

    logging.info('model {} created'.format(str(model)))

    np.random.seed(23456)
    cap = args.hop_cap
    transition = list(np.random.randn(2 * 2))

    nn_idx, efeature = generate_knn_table(args.chain_length, args.neighbour)

    if args.use_cuda:
        nn_idx, efeature = nn_idx.cuda(), efeature.cuda()
        model.cuda()
        emodel.cuda()

    # train_data_set = lib.data.RandomPGM(args.chain_length, cap, transition)
    # dataloader = torch.utils.data.DataLoader(train_data_set,
    #                                          batch_size=args.batch_size,
    #                                          shuffle=True,
    #                                          num_workers=8,
    #                                          worker_init_fn=worker_init_fn)

    train_dataset = lib.data.RandomPGMData(
        args.train_path, pgm_type="raw", size=args.train_size)
    test_dataset = lib.data.RandomPGMData(
        args.test_path, pgm_type="raw", size=args.test_size)

    train_loader = torch.utils.data.DataLoader(train_dataset,
                                               batch_size=args.batch_size,
                                               shuffle=True,
                                               num_workers=8,
                                               worker_init_fn=worker_init_fn)
    test_loader = torch.utils.data.DataLoader(test_dataset,
                                              batch_size=args.batch_size,
                                              shuffle=False,
                                              num_workers=8,
                                              worker_init_fn=worker_init_fn)

    optimizer = torch.optim.Adam(list(model.parameters()) +
                                 list(emodel.parameters()),
                                 lr=3e-3)
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_lambda=lambda x: max(0.98**x, 1e-6))
    start_epoch = 0
    gcnt = 0
    if os.path.exists(args.model_path):
        ckpt = torch.load(args.model_path)
        model.load_state_dict(ckpt['model_state_dict'])
        emodel.load_state_dict(ckpt['emodel_state_dict'])
        optimizer.load_state_dict(ckpt['optimizer_state_dict'])
        scheduler.load_state_dict(ckpt['lr_sche'])

        start_epoch = ckpt['epoch']
        gcnt = ckpt['gcnt']

    def get_model_dict():
        return {
                'model_state_dict': model.state_dict(),
                'emodel_state_dict': emodel.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'lr_sche': scheduler.state_dict(),
                'epoch': epoch,
                'gcnt': gcnt
            }

    def get_filename(epoch):
        return f'raw_nn_{args.neighbour}_epoches_{epoch}.pt'

    epoch = 0
    for epoch in tqdm(range(start_epoch, args.train_epoches)):
        torch.save(get_model_dict(), get_filename(epoch))

        logging.info(f'save train result to {get_filename(epoch)}')
        scheduler.step()

        loss_seq = []
        acc_seq = []
        acc_lp_seq = []
        for bcnt, (nfeature, nlabel, lp_label) in tqdm(enumerate(train_loader)):
            optimizer.zero_grad()
            if args.use_cuda:
                nfeature, nlabel, lp_label = nfeature.cuda(), nlabel.cuda(
                ), lp_label.cuda()
            if len(nfeature.shape) == 3:
                nfeature = nfeature.unsqueeze(-1)

            etype = emodel(efeature)
            # print(etype.shape)
            # print(nn_idx.shape)
            pred = model(nfeature, nn_idx.repeat(nfeature.shape[0], 1, 1),
                         etype.repeat(nfeature.shape[0], 1, 1, 1))
            pred = pred.squeeze(-1).permute(0, 2, 1).contiguous()
            loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                     nlabel.view(-1))
            loss.backward()
            torch.nn.utils.clip_grad_norm(
                list(model.parameters()) + list(emodel.parameters()), 1.0)

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
                writer.add_scalar('syn_train/loss', loss.item(), gcnt)
                writer.add_scalar('syn_train/acc', acc, gcnt)
                writer.add_scalar('syn_train/lp_acc', lp_acc, gcnt)
                logging.info(
                    'epoch = {} bcnt = {} loss = {} acc = {} lp_acc={}'.format(
                        epoch, bcnt, np.mean(loss_seq), np.mean(acc_seq),
                        np.mean(acc_lp_seq)))
                loss_seq = []
                acc_seq = []
                acc_lp_seq = []

    if epoch == args.train_epoches - 1:
        epoch = args.train_epoches
        torch.save(get_model_dict(), get_filename(epoch))

        logging.info(f'save train result to {get_filename(epoch)}')
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
    emodel.eval()
    for bcnt, (nfeature, nlabel,
                lp_label) in tqdm(enumerate(test_loader)):
        if args.use_cuda:
            nfeature, nlabel, lp_label \
                = nfeature.cuda(), nlabel.cuda(), lp_label.cuda()
        if len(nfeature.shape) == 3:
            nfeature = nfeature.unsqueeze(-1)

            etype = emodel(efeature)
            pred = model(nfeature, nn_idx.repeat(nfeature.shape[0], 1, 1),
                         etype.repeat(nfeature.shape[0], 1, 1, 1))
            pred = pred.squeeze(-1).permute(0, 2, 1).contiguous()
            loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                     nlabel.view(-1))
            torch.nn.utils.clip_grad_norm(
                list(model.parameters()) + list(emodel.parameters()), 1.0)

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

    logging.info(f'testing result: acc = {accum_acc / gcnt}, acc_lp = {accum_acc_lp / gcnt}')
    logging.info(f'stddev = {st.stdev(acc_global)}, stddev_lp = {st.stdev(acc_lp_global)}')

if __name__ == '__main__':
    main()
