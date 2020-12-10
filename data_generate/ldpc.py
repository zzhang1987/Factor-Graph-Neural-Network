import argparse
import numpy as np
import os
import torch
import lib.data.ldpc as ldpc
import lib.data.MNC as MNC
import logging
from multiprocessing import Pool



def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument('--num', type=int, default=1000)
    args.add_argument('--seed', type=int, default=42)
    args.add_argument('--save_path', type=str, required=True)
    args.add_argument('--burst_prob', type=float, default=0.05)
    return args.parse_args()

reseed = False

class ldpc_generate_worker:
    def __init__(self, snr_db, sigma_b, burst_prob):
        self.snr_db = snr_db
        self.sigma_b = sigma_b
        self.burst_prob = burst_prob

    def __call__(self, idx):
        global reseed
        if not reseed:
            np.random.seed(os.getpid())
            MNC.init_seed(os.getpid())
            logging.info('reseed at process: {}'.format(os.getpid()))
            reseed = True

        node_feature, hop_feature, gt, error_sum_product = \
            ldpc.gen_data_item(self.snr_db, self.sigma_b, self.burst_prob)
        return node_feature, hop_feature, gt, error_sum_product 
            

def generate_dataset():
    arg = parse_args()

    sigma_b = [0, 1, 2, 3, 4, 5]
    snr_db = [0, 1, 2, 3, 4]

    logging.basicConfig(filename='gen_ldpc_data.log', level=logging.INFO)
    error = np.zeros([len(snr_db), len(sigma_b)])

    node_features = []
    hop_features = []
    gts = []
    sigma_bs = []
    
    for j, csigma_b in enumerate(sigma_b):
        for i, csnr_db in enumerate(snr_db):
            logging.info(f'sigma_b = {csigma_b}, snr_db = {csnr_db}')

            worker = ldpc_generate_worker(csnr_db, csigma_b, arg.burst_prob)
            with Pool(24) as pool:
                data = pool.map(worker, list(range(arg.num)))

            node_feature = np.asarray([ditem[0] for ditem in data])
            hop_feature = np.asarray([ditem[1] for ditem in data])
            gt = np.asarray([ditem[2] for ditem in data], dtype=np.int)
            error_rate = np.mean(np.asarray([ditem[3] for ditem in data]))
            error[i, j] = error_rate

            node_features.append(node_feature)
            hop_features.append(hop_feature)
            gts.append(gt)
            sigma_bs.append(np.ones(arg.num) * csigma_b)

    node_features = torch.FloatTensor(np.concatenate(node_features, axis=0))
    hop_features = torch.FloatTensor(np.concatenate(hop_features, axis=0))
    sigma_bs = torch.FloatTensor(np.concatenate(sigma_bs, axis=0))
    gts = torch.LongTensor(np.concatenate(gts, axis=0))
    logging.info('LP Relaxation error on dataset = {}, avg_error = {}'.format(error, np.mean(error)))
    torch.save({
        'node_feature': node_features,
        'hop_feature': hop_features,
        'y': gts,
        'sigma_b': sigma_bs,
    }, arg.save_path)
    

    

            

            
    

if __name__ == '__main__':
    generate_dataset()
