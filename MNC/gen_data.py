import argparse
import numpy as np
import os
import random
import scipy.stats as st
import torch
from tqdm import tqdm


def get_arguments():
    args = argparse.ArgumentParser()
    args.add_argument('--train', action='store_true', default=False)
    args.add_argument('--num', type=int, default=100)
    args.add_argument('--seed', type=int, default=42)
    return args.parse_args()

def get_snr(snr_db):
    return 10**(snr_db/10)


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


def work(gcx, sigma, train):
    print(f'gcx = {gcx}, sigma = {sigma}')

    # 1. original code: s
    s = np.random.randint(0, 2, 48)
    write_to_file(s, '_s/48')

    # 2. transmitted code: t (groundtruth)
    os.system('./s2t -sfile _s/48 -k 48 -n 48 -Gfile codes/96.3.963/G -smn 1 -tfile _t/96 ')
    t = read_from_file('_t/96', int)

    # 3. received code: y
    os.system(
        f'./t2y -tfile _t/96 -yfile _y/96 -gcx {gcx} -seed 322457 -n 96 -sigma {sigma} ')
    y = read_from_file('_y/96', float)

    error = 0

    if not train:
        # 4. sum-product decoding
        os.system(f'./y2b -yfile _y/96 -n 96 -bfile _b/96 -gcx {gcx} ')
        os.system('./zb2x -bfile _b/96 -zfixed 0 -k 48 -n 48 -Afile codes/96.3.963/A2 -xfile _x/48gc -xso 1 -bndloops 100 ')
        x = read_from_file('_x/48gc', int)

        # 5. calculate BER
        for a, b in zip(s, x):
            if a != b:
                error += 1

    return y, np.ones((9, 48)), t, error


def gen_value():
    args = get_arguments()

    sigma_b = [0, 1, 2, 3, 4, 5]
    snr_db = [0, 1, 2, 3, 4]

    snr = [get_snr(x) for x in snr_db]
    sigma_c = [1 / x for x in snr]

    l = []

    np.random.seed(args.seed)

    for i, x in enumerate(sigma_b):
        for j, y in enumerate(sigma_c):
            sigma_a = x / y
            print(f'sigma_b = {x}, snr_db = {snr_db[j]}, snr = {snr[j]}, sigma_a = {sigma_a}')
            l.append((x, snr_db[j], snr[j], sigma_a))

    node_feature = []
    hop_feature = []
    y_list = []
    sigma_b_list = []

    if not args.train:
        error = [[[] for j in range(6)] for i in range(5)]
        err_mean = np.zeros((5, 6))
        err_err = np.zeros((5, 6))

    for i in tqdm(range(5)):
        for _ in tqdm(range(args.num)):
            j = random.randint(0, 5)
            idx = j * 5 + i
            x, m, y, err = work(l[idx][2], l[idx][3], args.train)

            node_feature.append([[elem, l[idx][2]] for elem in x])
            hop_feature.append(m)
            y_list.append(y)
            sigma_b_list.append(j)

            if not args.train:
                error[i][j].append(err)

    node_feature = torch.FloatTensor(np.stack(node_feature)).permute(0, 2, 1).unsqueeze(-1)
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

    error /= 48

    if not args.train:
        for i in range(5):
            for j in range(6):
                err_mean[i][j] = np.mean(error[i][j])
                err_err[i][j] = st.sem(error[i][j])

        print(err_mean)
        print(err_err)

        torch.save({
            'err_mean': err_mean,
            'err_err': err_err,
        }, 'err.pt')


gen_value()
