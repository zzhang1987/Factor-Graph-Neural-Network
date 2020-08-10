import torch
import numpy as np
import time


def generate_knn_table(n: int,
                       k: int,
                       generate_efeature: bool = False):
    """generate graph structure for a chain model

    :param n: number of nodes
    :param k: window size
    :param generate_efeature: generate edge feature or not
    :returns: 
    :rtype: 

    """
    if k % 2 == 0:
        k = k + 1
    nn_idx = np.zeros([n, k]).astype(np.int64)
    if generate_efeature:
        efeature = np.zeros([1, n, k]).astype(np.float32)
    hk = k // 2
    for i in range(n):
        for idx, j in enumerate(range(i - hk, i + hk)):
            if j < 0:
                j = 0
            if j >= n:
                j = n - 1
            nn_idx[i, idx] = j
            if generate_efeature:
                efeature[0, i, idx] = i - j
    nn_idx = torch.from_numpy(np.expand_dims(nn_idx, 0))
    if generate_efeature:
        efeature = torch.from_numpy(np.expand_dims(efeature, 0))
        return nn_idx, efeature
    else:
        return nn_idx


def generate_pw_factor_table(n):
    """generate pairwise graph structrure for a chain

    :param n: number of nodes
    :returns: 
    :rtype: 

    """
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
    """generate higher order graph structure

    :param n: number of nodes
    :param k: window size
    :returns: 
    :rtype: 

    """
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
