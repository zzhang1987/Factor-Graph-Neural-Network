import model.syn_models as syn_models
from utils.syn_utils import generate_high_factor_table, generate_knn_table, generate_pw_factor_table, worker_init_fn
import lib
import torch
import argparse
import sys
import utils
import logging
import datetime
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
    parser.add_argument('--legacy_model_path', type=str)
    parser.add_argument('--model_name',
                        type=str,
                        default='CHOPModel',
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


def main():
    args = parse_args()
    subdir = f'train_syn_hop_factor_{args.model_name}_nn_{args.neighbour}_at_{datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S")}'
    logdir = f'./tf_logs/{subdir}'
    writer = SummaryWriter(log_dir=logdir)
    utils.init_logger('./logs/', subdir, print_log=True)
    logging.info(str(args))

    nfeature_dim = 2
    model_create_eval = 'syn_models.{}({}, {})'.format(
        args.model_name, nfeature_dim, args.hop_order)
    model = eval(model_create_eval)
    logging.info('model {} created'.format(model))

    if hasattr(args, 'legacy_model_path'):
        if args.legacy_model_path is not None and os.path.exists(args.legacy_model_path):
            model.load_legacy_weight(args.legacy_model_path)

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
    train_dataset = lib.data.RandomPGMData(args.train_path,
                                           pgm_type="hops",
                                           size=args.train_size)
    test_dataset = lib.data.RandomPGMData(args.test_path,
                                          pgm_type="hops",
                                          size=args.test_size)

    train_loader = torch.utils.data.DataLoader(train_dataset,
                                               batch_size=args.batch_size,
                                               shuffle=True,
                                               num_workers=16,
                                               worker_init_fn=worker_init_fn)
    test_loader = torch.utils.data.DataLoader(test_dataset,
                                              batch_size=args.batch_size,
                                              shuffle=True,
                                              num_workers=16,
                                              worker_init_fn=worker_init_fn)

    optimizer = torch.optim.Adam(model.parameters(), lr=3e-3)
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer, lr_lambda=lambda x: max(0.998**x, 1e-6))
    start_epoch = 0
    gcnt = 0

    if os.path.exists(args.model_path):
        ckpt = torch.load(args.model_path)
        model.load_state_dict(ckpt['model_state_dict'])
        optimizer.load_state_dict(ckpt['optimizer_state_dict'])
        scheduler.load_state_dict(ckpt['lr_sche'])
        start_epoch = ckpt['epoch']
        gcnt = ckpt['gcnt']
        print('model loaded')

    def get_model_dict():
        return {
            'model_state_dict': model.state_dict(),
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

            pred = model(nfeature, pws, hops, nn_idx_pw, nn_idx_high,
                         efeature_pw, efeature_high)
            loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                     nlabel.view(-1))
            loss.backward()
            torch.nn.utils.clip_grad_norm(model.parameters(), 1.0)
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
    for bcnt, (nfeature, pws, hops, nlabel,
               lp_label) in tqdm(enumerate(test_loader)):
        if args.use_cuda:
            nfeature, pws, hops, nlabel, lp_label \
                = nfeature.cuda(), pws.cuda(), hops.cuda(), nlabel.cuda(), lp_label.cuda()
        if len(nfeature.shape) == 3:
            nfeature = nfeature.unsqueeze(-1)

        with torch.no_grad():
            pred = model(nfeature, pws, hops, nn_idx_pw, nn_idx_high,
                         efeature_pw, efeature_high)
            loss = torch.nn.functional.cross_entropy(pred.view(-1, 2),
                                                     nlabel.view(-1))

        loss_seq.append(loss.item())
        gcnt += 1
        # print(pred.shape)
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
