import argparse
import datetime
import logging
import numpy as np
import os
import sys
import torch
from tensorboardX import SummaryWriter
import time
from tqdm import tqdm

import lib
from lib.mpnn import factor_mpnn
import scipy.stats as st
import utils
from utils.types import str2bool


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--n_epochs',
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
    parser.add_argument('--use_cuda',
                        type=str2bool,
                        default=True,
                        help="Use cuda or not")

    parser.add_argument('--train_path',
                        type=str,
                        default="ldpc_data/train.pt",
                        help="path of the training dataset")

    parser.add_argument('--test_path',
                        type=str,
                        default="ldpc_data/valid.pt",
                        help="path of the testing dataset")

    parser.add_argument('--filename', type=str, default='ldpc_data/96.3.963')
    parser.add_argument('--train', action='store_true', default=False)

    parser.add_argument('--batch_size', type=int, default=32)

    return parser.parse_args()


def generate_high_factor_table(filename):
    n_nodes = 96+48
    n_edges = 6+3
    nn_idx = np.zeros((n_nodes, n_edges)).astype(np.int64)
    efeature = np.zeros((2, n_nodes, n_edges)).astype(np.float32)

    with open(filename) as f:
        for i in range(4):
            f.readline()
        for i in range(96):
            indice = map(lambda x: int(x)-1, f.readline().strip().split())
            for j, idx in enumerate(indice):
                nn_idx[i, j] = 96 + idx
                efeature[0, i, j] = 1

        for i in range(48):
            indice = map(lambda x: int(x)-1, f.readline().strip().split())
            for j, idx in enumerate(indice):
                nn_idx[96+i, 3+j] = idx
                efeature[1, 96+i, 3+j] = 1

    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    efeature = torch.from_numpy(np.expand_dims(efeature, 0))
    return nn_idx, efeature


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    np.random.seed(((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) +
                   ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24))


def train(args, model, emodel_high, nn_idx_high, efeature_high, writer, model_dir):
    #train_dataset = lib.data.Codes(args.train_path)

    train_dataset = lib.data.ContinusCodes()

    train_loader = torch.utils.data.DataLoader(train_dataset,
                                               batch_size=args.batch_size,
                                               shuffle=True,
                                               num_workers=8,
                                               worker_init_fn=worker_init_fn)

    parameters = list(model.parameters()) + list(emodel_high.parameters())
    optimizer = torch.optim.Adam(parameters, lr=1e-4, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_lambda=lambda x: max(0.98**x, 1e-6))
    start_epoch = 0
    gcnt = 0
    if os.path.exists(args.model_path):
        ckpt = torch.load(args.model_path)
        model.load_state_dict(ckpt['model_state_dict'])
        emodel_high.load_state_dict(ckpt['emodel_high_state_dict'])
        optimizer.load_state_dict(ckpt['optimizer_state_dict'])
        scheduler.load_state_dict(ckpt['lr_sche'])

        start_epoch = ckpt['epoch']
        gcnt = ckpt['gcnt']

    def get_model_dict():
        return {
            'model_state_dict': model.state_dict(),
            'emodel_high_state_dict': emodel_high.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'lr_sche': scheduler.state_dict(),
            'epoch': epoch,
            'gcnt': gcnt
        }

    def get_model_path():
        return os.path.join(model_dir, '{}_nn_factor_epoches_{}.pt'.format(args.model_name, epoch))

    print('training started!')
    for epoch in tqdm(range(start_epoch, args.n_epochs)):
        torch.save(get_model_dict(), get_model_path())

        logging.info(f'save train result to {get_model_path()}')

        loss_seq = []
        acc_seq = []
        for bcnt, (nfeature, hops, label, _) in tqdm(enumerate(train_loader)):
            optimizer.zero_grad()
            if args.use_cuda:
                # print(nfeature)
                # print(hops)
                # print(label)
                nfeature, hops, label \
                    = nfeature.cuda(), hops.cuda(), label.cuda()
            hops = hops.float()
            if len(nfeature.shape) == 3:
                nfeature = nfeature.unsqueeze(-1)

            etype_high = emodel_high(efeature_high)
            bsize = nfeature.shape[0]

            # print('efeature_high', efeature_high.shape)
            # print('etype_high', etype_high.shape)

            # print('nfeature', nfeature.shape)
            # print('hops', hops.shape)

            # print('nn_idx_high', nn_idx_high.shape)
            # print('etype_high', etype_high.shape)

            pred, _ = model(nfeature, [hops],
                            [[
                                nn_idx_high.repeat(bsize, 1, 1),
                                etype_high.repeat(bsize, 1, 1, 1)
                            ]])
            # print(pred.shape)
            # print(label.shape)

            pred = pred.squeeze()[:, :48].contiguous()
            label = label[:, :48].contiguous()

            # print(label.shape)
            loss = torch.nn.functional.binary_cross_entropy_with_logits(
                pred.view(-1), label.view(-1).float())
            loss.backward()
            # torch.nn.utils.clip_grad_norm(parameters, 1.0)

            optimizer.step()
            loss_seq.append(loss.item())
            gcnt += 1

            pred_int = (pred[:, :48] > 0)
            all_correct = torch.sum(pred_int.long() == label)
            acc = all_correct.item() / np.prod(label.shape)

            acc_seq.append(acc)

            if gcnt % 10 == 0:
                logging.info('epoch = {} bcnt = {} loss = {} acc = {}'.format(
                    epoch, bcnt, np.mean(loss_seq), np.mean(acc_seq)))
                writer.add_scalar('syn_train/loss', np.mean(loss_seq), gcnt)
                writer.add_scalar('syn_train/acc', np.mean(acc_seq), gcnt)
                loss_seq = []
                acc_seq = []

        scheduler.step()

    if epoch == args.n_epochs - 1:
        epoch = args.n_epochs
        torch.save(get_model_dict(), get_model_path())
        logging.info(f'save train result to {get_model_path()}')
        logging.info('training done!')


def test(args, model, emodel_high, nn_idx_high, efeature_high):
    test_dataset = lib.data.Codes(args.test_path, train=False)

    test_loader = torch.utils.data.DataLoader(test_dataset,
                                              batch_size=10,
                                              shuffle=False,
                                              num_workers=8,
                                              worker_init_fn=worker_init_fn)

    assert args.model_path, 'please input model path'

    ckpt = torch.load(args.model_path)
    model.load_state_dict(ckpt['model_state_dict'])
    emodel_high.load_state_dict(ckpt['emodel_high_state_dict'])

    acc_seq = []
    acc_cnt = np.zeros((5, 6))
    acc_tot = np.zeros((5, 6))
    tot = 0
    model.eval()
    emodel_high.eval()

    SNR = [0, 1, 2, 3, 4]
    for _, (nfeature, hops, label, sigma_b) in tqdm(enumerate(test_loader)):
        if args.use_cuda:
            # print(nfeature)
            # print(hops)
            # print(label)
            # print(sigma_b)
            nfeature, hops, label, sigma_b \
                = nfeature.cuda(), hops.cuda(), label.cuda(), sigma_b.cuda()
        cur_SNR = nfeature[:, 1, 0, 0]
        # print(cur_SNR)
        hops = hops.float()
        # print(cur_SNR)

        if len(nfeature.shape) == 3:
            nfeature = nfeature.unsqueeze(-1)

        etype_high = emodel_high(efeature_high)
        bsize = nfeature.shape[0]

        pred, _ = model(
            nfeature, [hops],
            [[
                nn_idx_high.repeat(bsize, 1, 1),
                etype_high.repeat(bsize, 1, 1, 1)
            ]])

        pred = pred.squeeze().contiguous()
        pred_int = (pred > 0).long()

        for i, elem in enumerate(SNR):
            for b in range(6):
                indice = (sigma_b == b) & (abs(cur_SNR-elem) < 1e-3)
                acc_cnt[i][b] += torch.sum(pred_int[indice, :48]
                                           == label[indice, :48])
                acc_tot[i][b] += torch.sum(indice) * 48

        # parameters = list(model.parameters()) + list(emodel_high.parameters())
        # torch.nn.utils.clip_grad_norm(parameters, 1.0)

        all_correct = torch.sum(pred_int[:, :48] == label[:, :48])

        indice = sigma_b
        acc_seq.append(all_correct.item())
        tot += np.prod(label.shape) // 2

    print(1 - sum(acc_seq) / tot)
    err_class = 1 - np.divide(acc_cnt, acc_tot)
    print(torch.FloatTensor(err_class))


def main():
    args = parse_args()

    nfeature_dim = 2
    hop_order = 6
    if args.model_name == 'mp_nn_factor':
        model = factor_mpnn(nfeature_dim, [hop_order],
                            [64, 64, 128, 128, 256, 256, 128, 128, 64, 64, 1],
                            [8])

        emodel_high = torch.nn.Sequential(torch.nn.Conv2d(2, 64, 1),
                                          torch.nn.ReLU(inplace=True),
                                          torch.nn.Conv2d(64, 8, 1))

    def get_model_description():
        return str(model) + str(emodel_high)

    logging.info('model {} created'.format(get_model_description()))

    nn_idx_high, efeature_high = generate_high_factor_table(args.filename)

    if args.use_cuda:
        nn_idx_high = nn_idx_high.cuda()
        efeature_high = efeature_high.cuda()

        model.cuda()
        emodel_high.cuda()

    if args.train:
        subdir = f'train_syn_hop_factor_{args.model_name}_at_{datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}'
        utils.init_logger('./logs/', subdir, print_log=True)
        logging.info(str(args))
        logdir = f'./tf_logs/{subdir}'
        print(f'logdir = {logdir}')
        writer = SummaryWriter(log_dir=logdir)

        model_dir = f'./model_ldpc/{subdir}'
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)
        train(args, model, emodel_high, nn_idx_high,
              efeature_high, writer, model_dir)

    else:
        test(args, model, emodel_high, nn_idx_high, efeature_high)


if __name__ == '__main__':
    main()
