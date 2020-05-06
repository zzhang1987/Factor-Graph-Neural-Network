import argparse
import numpy as np
import os
import random
import scipy.stats as st
import torch
from tqdm import tqdm
from MNC import s2t, t2y


def get_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--train', action='store_true', default=False)
    args.add_argument('--num', type=int, default=100)
    args.add_argument('--seed', type=int, default=42)
    return args.parse_args()


def get_snr(snr_db):
    return 10**(snr_db/20)


def write_to_file(s, filename):
    with open(filename, 'w') as f:
        for elem in s:
            f.write(str(elem)+'\n')


def read_from_file(filename, type):
    ret = []
    with open(filename) as f:
        for line in f:
            ret.append(type(line.strip()))
    return ret


def work(snr_db, sigma_b, train):
    print(f'snr_db = {snr_db}, sigma_b = {sigma_b}')

    gcx = get_snr(snr_db)
    # 1. original code: s
    s = np.random.randint(0, 2, 48)
    # write_to_file(s, '_s/48')

    # # 2. transmitted code: t (groundtruth)
    # os.system(
    #     './s2t -sfile _s/48 -k 48 -n 48 -Gfile codes/96.3.963/G -smn 1 -tfile _t/96 ')
    # t = read_from_file('_t/96', int)
    t = s2t(s, 48, 48, './codes/96.3.963/G', True)

    y = y = t2y(t, snr_db, sigma_b, 0.05)

    write_to_file(y, '_y/96')
    error = error_prime = error_b = 0

    if not train:
        # 4.1 [normal LDPC] sum-product decoding
        os.system(f'./y2b -yfile _y/96 -n 96 -bfile _b/96 -gcx {gcx} ')
        os.system(
            './zb2x -bfile _b/96 -zfixed 0 -k 48 -n 48 -Afile codes/96.3.963/A2 -xfile _x/48gc -xso 1 -bndloops 100 ')
        x = read_from_file('_x/48gc', int)

        # 4.2 [normal LDPC] calculate BER
        error = sum(a != b for a, b in zip(s, x))

        # 5. [bursty LDPC] sum-product decoding
        gcx_prime = 1 / np.sqrt(1/(gcx**2) + 12/5)
        os.system(f'./y2b -yfile _y/96 -n 96 -bfile _b/96 -gcx {gcx_prime} ')
        os.system(
            './zb2x -bfile _b/96 -zfixed 0 -k 48 -n 48 -Afile codes/96.3.963/A2 -xfile _x/48gc -xso 1 -bndloops 100 ')
        x = read_from_file('_x/48gc', int)
        error_prime = sum(a != b for a, b in zip(s, x))

        # 6. [bits baseline]
        error_b = 0
        for a, b in zip(t, y):
            p1 = st.norm(gcx, 1).pdf(b)
            p2 = st.norm(-gcx, 1).pdf(b)
            error_b += (p1 >= p2) != a

    return y, np.ones((9, 48)), t, error, error_prime, error_b


def process_error(s, error):
    print(f'======       {s}        ======')
    err_mean = np.zeros((5, 6))
    err_err = np.zeros((5, 6))

    for i in range(5):
        for j in range(6):
            err_mean[i][j] = np.mean(error[i][j])
            err_err[i][j] = st.sem(error[i][j])

    err_mean /= 48
    err_err /= 48**2

    print(err_mean)
    print(err_err)

    torch.save({
        'err_mean': err_mean,
        'err_err': err_err,
    }, f'err_{s}.pt')


def gen_value():
    args = get_arguments()

    sigma_b = [0, 1, 2, 3, 4, 5]
    snr_db = [0, 1, 2, 3, 4]

    l = []

    np.random.seed(args.seed)
    node_feature = []
    hop_feature = []
    y_list = []
    sigma_b_list = []

    if not args.train:
        error = [[[] for j in range(6)] for i in range(5)]
        error_prime = [[[] for j in range(6)] for i in range(5)]
        error_b = [[[] for j in range(6)] for i in range(5)]

    for j, csigma_b in enumerate(sigma_b):
        for i, csnr_db in enumerate(snr_db):
            print(
                f'sigma_b = {csigma_b}, snr_db = {csnr_db}')

            for num in range(args.num):
                x, m, y, err, err_prime, err_b = work(
                    csnr_db, csigma_b, args.train)
                node_feature.append([[elem, csnr_db] for elem in x])
                hop_feature.append(m)
                y_list.append(y)
                sigma_b_list.append(j)
                if not args.train:
                    error[i][j].append(err)
                    error_prime[i][j].append(err_prime)
                    error_b[i][j].append(err_b)

    node_feature = torch.FloatTensor(
        np.stack(node_feature)).permute(0, 2, 1).unsqueeze(-1)
    hop_feature = torch.LongTensor(np.stack(hop_feature)).unsqueeze(-1)
    y = torch.LongTensor(np.stack(y_list))
    sigma_b = torch.LongTensor(np.stack(sigma_b_list))

    print(node_feature.shape)
    print(hop_feature.shape)
    print(y.shape)
    print(sigma_b.shape)

    filename = 'train.pt' if args.train else 'valid.pt'
    torch.save({
        'node_feature': node_feature,
        'hop_feature': hop_feature,
        'y': y,
        'sigma_b': sigma_b,
    }, filename)

    if not args.train:
        process_error('ldpc', error)
        process_error('ldpc-bursty', error_prime)
        process_error('bits-baseline', error_b)


gen_value()
