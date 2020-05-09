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
from lib.mpnn import factor_mpnn, FactorNN
import scipy.stats as st
import utils
from utils.types import str2bool, to_cuda


class LDPCModel(torch.nn.Module):
    def __init__(self, nfeature_dim, hop_order, nedge_type, with_residual=True):
        super(LDPCModel, self).__init__()

        self.main = FactorNN(nfeature_dim,
                             [hop_order, 96],
                             [64, 64, 64, 128, 128, 64, 64],
                             [nedge_type, 1],
                             2,
                             skip_link={3: 2, 4: 1, 5: 0},
                             ret_high=True)

        self.emodel_f2v = torch.nn.Sequential(torch.nn.Conv2d(7, 64, 1),
                                              torch.nn.ReLU(inplace=True),
                                              torch.nn.Conv2d(64, nedge_type, 1))

        self.emodel_v2f = torch.nn.Sequential(torch.nn.Conv2d(7, 64, 1),
                                              torch.nn.ReLU(inplace=True),
                                              torch.nn.Conv2d(64, nedge_type, 1))

        hetype_v2f = np.ones([1, 1, 1, 96]).astype(np.float32)
        hetype_f2v = np.ones([1, 1, 96, 1]).astype(np.float32)

        hnn_idx_v2f = np.asarray(
            list(range(96)), dtype=np.int64).reshape(1, 1, 96)
        hnn_idx_f2v = np.asarray([0] * 96, dtype=np.int64).reshape(1, 96, 1)

        self.hnn_idx_v2f = torch.nn.Parameter(
            torch.from_numpy(hnn_idx_v2f).long(), requires_grad=False)
        self.hnn_idx_f2v = torch.nn.Parameter(
            torch.from_numpy(hnn_idx_f2v).long(), requires_grad=False)

        self.hetype_v2f = torch.nn.Parameter(
            torch.from_numpy(hetype_v2f), requires_grad=False)
        self.hetype_f2v = torch.nn.Parameter(
            torch.from_numpy(hetype_f2v), requires_grad=False)

        self.with_residual = with_residual

        self.nhop_regressor = torch.nn.Sequential(torch.nn.Linear(64, 128),
                                                  torch.nn.BatchNorm1d(128),
                                                  torch.nn.ReLU(),
                                                  torch.nn.Linear(128, 128),
                                                  torch.nn.ReLU(),
                                                  torch.nn.Linear(128, 1),
                                                  torch.nn.ReLU())

    def forward(self, node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f):
        etype_f2v = self.emodel_f2v(efeature_f2v)
        etype_v2f = self.emodel_v2f(efeature_v2f)

        with torch.no_grad():
            bsize = node_feature.shape[0]
            nhop_feature = node_feature[:, 0, :, :]
            nhop_feature = nhop_feature.reshape(bsize, 96, 1, 1)

        # print(nn_idx_f2v[0, :, :].shape)
        # print(nn_idx_v2f[0, :, :].shape)
        res, nhops = self.main(node_feature,
                               [hop_feature, nhop_feature],
                               [nn_idx_f2v, self.hnn_idx_f2v.repeat(
                                   bsize, 1, 1)],
                               [nn_idx_v2f, self.hnn_idx_v2f.repeat(
                                   bsize, 1, 1)],
                               [etype_f2v, self.hetype_f2v.repeat(
                                   bsize, 1, 1, 1)],
                               [etype_v2f, self.hetype_v2f.repeat(bsize, 1, 1, 1)])

        if self.with_residual:
            res = res + node_feature[:, :1, :, :]
        batch_size = res.shape[0]
        res = res.squeeze()
        if batch_size == 1:
            res = res.unsqueeze(0)

        hhop = nhops[1].squeeze()
        if bsize == 1:
            hhop = hhop.unsqueeze(0)

        snr_b_pred = self.nhop_regressor(hhop)

        return res[:, :48].contiguous(), snr_b_pred


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
                        default='FactorNN',
                        help="model name")
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
    parser.add_argument('--train', action='store_true', default=True)

    parser.add_argument('--batch_size', type=int, default=32)
    parser.add_argument('--aggregator', type=str, default='max')
    args = parser.parse_args()

    if not torch.cuda.is_available():
        args.use_cuda = False

    return args


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    seed = ((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) + \
        ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24)
    np.random.seed(seed)
    lib.data.init_seed(int(seed % 65537))


def train(args, model,  writer, model_dir):

    train_dataset = lib.data.ContinousCodesSP()

    train_loader = torch.utils.data.DataLoader(train_dataset,
                                               batch_size=args.batch_size,
                                               shuffle=True,
                                               num_workers=8,
                                               worker_init_fn=worker_init_fn)

    optimizer = torch.optim.Adam(
        model.parameters(), lr=1e-2,  weight_decay=1e-8)

    def lr_sched(x, start=10):
        if x <= start:
            return max(1e-2, (1.0 / start) * x)
        else:
            return max(0.99 ** (x - start), 1e-6)
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_lambda=lambda x: lr_sched(x))
    start_epoch = 0
    gcnt = 0
    if os.path.exists(args.model_path):
        if not torch.cuda.is_available():
            ckpt = torch.load(
                args.model_path, map_location=torch.device('cpu'))
        else:
            ckpt = torch.load(args.model_path)
        model.load_state_dict(ckpt['model_state_dict'])
        optimizer.load_state_dict(ckpt['optimizer_state_dict'])
        scheduler.load_state_dict(ckpt['lr_sche'])

        start_epoch = ckpt['epoch']
        gcnt = ckpt['gcnt']

    def get_model_dict():
        return {
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'lr_sche': scheduler.state_dict(),
            'epoch': epoch,
            'gcnt': gcnt
        }

    def get_model_path():
        return os.path.join(model_dir, '{}_nn_factor_epoches_{}.pt'.format(args.model_name, epoch))

    print('training started!')
    for epoch in tqdm(range(start_epoch, args.n_epochs)):
        if (epoch + 1) % 10 == 0:
            torch.save(get_model_dict(), get_model_path())

        logging.info(f'save train result to {get_model_path()}')

        loss_seq = []
        sigma_b_loss_seq = []
        acc_seq = []
        for bcnt, (node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label, sigma_b) in tqdm(enumerate(train_loader)):
            optimizer.zero_grad()
            if args.use_cuda:
                node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label, sigma_b = to_cuda(
                    node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label, sigma_b.float())

            if len(node_feature.shape) == 3:
                node_feature = node_feature.unsqueeze(-1)

            pred, sigma_b_pred = model(node_feature, hop_feature, nn_idx_f2v,
                                       nn_idx_v2f, efeature_f2v, efeature_v2f)

            label = label[:, :48].contiguous()
            # print(pred.shape)
            # print(label.shape)
            loss = torch.nn.functional.binary_cross_entropy_with_logits(
                pred.view(-1), label.view(-1).float())
            sigma_b_loss = torch.nn.functional.mse_loss(
                sigma_b_pred.view(-1), (torch.pow(10.0, sigma_b.float()/20)).view(-1))

            allloss = loss + 0.1 * sigma_b_loss
            allloss.backward()
            # torch.nn.utils.clip_grad_norm(parameters, 1.0)

            optimizer.step()
            loss_seq.append(loss.item())
            gcnt += 1

            pred_int = (pred > 0)
            all_correct = torch.sum(pred_int.long() == label)
            acc = all_correct.item() / np.prod(label.shape)

            acc_seq.append(acc)
            sigma_b_loss_seq.append(sigma_b_loss.item())

            if gcnt % 10 == 0:
                logging.info('epoch = {} bcnt = {} loss = {} acc = {}'.format(
                    epoch, bcnt, np.mean(loss_seq), np.mean(acc_seq)))
                writer.add_scalar('syn_train/loss', np.mean(loss_seq), gcnt)
                writer.add_scalar('syn_train/sigma_b_loss',
                                  np.mean(sigma_b_loss_seq), gcnt)
                writer.add_scalar('syn_train/acc', np.mean(acc_seq), gcnt)
                loss_seq = []
                acc_seq = []
                sigma_b_loss_seq = []

        scheduler.step()

    if epoch == args.n_epochs - 1:
        epoch = args.n_epochs
        torch.save(get_model_dict(), get_model_path())
        logging.info(f'save train result to {get_model_path()}')
        logging.info('training done!')


def test(args, model):
    test_dataset = lib.data.Codes_SP(args.test_path, train=False)

    test_loader = torch.utils.data.DataLoader(test_dataset,
                                              batch_size=1,
                                              shuffle=False,
                                              num_workers=8,
                                              worker_init_fn=worker_init_fn)

    assert args.model_path, 'please input model path'

    if args.use_cuda:
        ckpt = torch.load(args.model_path)
    else:
        ckpt = torch.load(args.model_path, map_location=torch.device('cpu'))

    model.load_state_dict(ckpt['model_state_dict'])

    acc_seq = []
    acc_cnt = np.zeros((5, 6))
    acc_tot = np.zeros((5, 6))
    tot = 0
    model.eval()
    # pdb.set_trace()
    SNR = [0, 1, 2, 3, 4]
    for _, (node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label, sigma_b) in enumerate(test_loader):
        if args.use_cuda:
            # print(nfeature)
            # print(hops)
            # print(label)
            # print(sigma_b)
            if args.use_cuda:
                node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label = to_cuda(
                    node_feature, hop_feature, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f, label)

        if len(node_feature.shape) == 3:
            node_feature = node_feature.unsqueeze(-1)
        cur_SNR = node_feature[:, 1, 0, 0]

        with torch.no_grad():
            pred, _ = model(node_feature, hop_feature, nn_idx_f2v,
                            nn_idx_v2f, efeature_f2v, efeature_v2f)

            pred = pred.squeeze().contiguous()
        pred_int = (pred >= 0).long().squeeze()
        label = label.squeeze()

        b = sigma_b.squeeze().item()
        i = round(cur_SNR.squeeze().item())

        # print(i, b)
        acc_cnt[i][b] += torch.sum(pred_int[:48] == label[:48]).item()
        acc_tot[i][b] += 48
        print(
            'snr = {} sigma_b = {}  Correct = {}/48'.format(i, b, torch.sum(pred_int[:48] == label[:48]).item()))

        # parameters = list(model.parameters()) + list(emodel_high.parameters())
        # torch.nn.utils.clip_grad_norm(parameters, 1.0)

        all_correct = torch.sum(pred_int[:48] == label[:48])

        acc_seq.append(all_correct.item())
        tot += np.prod(label.shape) // 2

    print(1 - sum(acc_seq) / tot)
    err_class = 1 - np.divide(acc_cnt, acc_tot)
    print(torch.FloatTensor(err_class))


def residual_link(final_feature, original_feature):
    return final_feature + original_feature[:, :1, :, :]


def main():
    args = parse_args()

    nfeature_dim = 2
    hop_order = 6
    nedge_types = 8
    model = LDPCModel(nfeature_dim, hop_order, nedge_types)

    def get_model_description():
        return str(model)

    logging.info('model {} created'.format(get_model_description()))

    if args.use_cuda:

        model.cuda()

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
        train(args, model, writer, model_dir)

    else:
        test(args, model)


if __name__ == '__main__':
    main()
