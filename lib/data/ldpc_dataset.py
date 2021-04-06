import torch
from torch.utils.data import Dataset

import numpy as np
import os

from .ldpc import gen_data_item
from .MNC import init_seed


class ldpc_graph_structure_generator:
    def __init__(self):
        self.code_len = 48
        self.ldpc_file = os.path.join(
            os.path.dirname(__file__),
            '../../ldpc_codes/96.3.963/96.3.963'
        )
        self.get_factor_structure()

    def get_factor_structure(self):
        n_nodes = 96 + 48
        n_edges = 6
        nn_idx = np.zeros((n_nodes, n_edges)).astype(np.int64)
        efeature = np.zeros((2, n_nodes, n_edges)).astype(np.float32)

        with open(self.ldpc_file) as f:
            for i in range(4):
                f.readline()

            for i in range(96):
                indice = map(lambda x: int(x)-1, f.readline().strip().split())
                for j, idx in enumerate(indice):
                    nn_idx[i, j] = 96 + idx
                    efeature[0, i, j] = 1

                for k in range(j + 1, n_edges):
                    nn_idx[i, k] = i
            # print(nn_idx[:96, :])
            factors = []

            for i in range(48):
                indice = map(lambda x: int(x)-1, f.readline().strip().split())
                cfactor = []
                for j, idx in enumerate(indice):
                    cfactor.append(idx)
                    nn_idx[96 + i, j] = idx
                    efeature[1, 96 + i, j] = 1

                factors.append(cfactor)

            self.factors = np.asarray(factors, np.int)
            self.etype = efeature
            self.nn_idx = nn_idx

    def get_pairwise_structure(self):
        n_nodes = 96
        n_edges = 3 * 6  # upper bound
        nn_idx = np.zeros((n_nodes, n_edges)).astype(np.int64)
        efeature = np.zeros((1, n_nodes, n_edges)).astype(np.float32)

        v2f = np.zeros((96, 3)).astype(np.int64)
        f2v = np.zeros((48, 6)).astype(np.int64)

        with open(self.ldpc_file) as f:
            for i in range(4):
                f.readline()

            for i in range(96):
                v2f[i, :] = list(
                    map(lambda x: int(x)-1, f.readline().strip().split()))

            for i in range(48):
                f2v[i, :] = list(
                    map(lambda x: int(x)-1, f.readline().strip().split()))

            for i in range(96):
                k = 0
                for f in v2f[i]:
                    for v in f2v[f]:
                        nn_idx[i, k] = v
                        efeature[0, i, k] = 1
                        k += 1

            self.factors = f2v
            self.etype = efeature
            self.nn_idx = nn_idx

    def get_highorder_feature(self, y):

        return np.take(y, self.factors.reshape(-1)).reshape(self.factors.shape)

    def get_mpnn_sp_structure(self, y):
        hop = self.get_highorder_feature(y)
        # hop = torch.sum(hop, dim=1)
        nn_idx_f2v = self.nn_idx[:96, :3] - 96
        nn_idx_v2f = self.nn_idx[96:, :]
        efeature_f2v = np.take(hop, nn_idx_f2v.reshape(-1),
                               axis=0).reshape(96, 3, 6).astype(np.float32)
        efeature_f2v = np.concatenate(
            [efeature_f2v, np.repeat(y.reshape(96, 1, 1), 3, axis=1)], axis=2)

        efeature_v2f = np.repeat(hop.reshape(48, 1, 6), 6, axis=1)
        efeature_v2f = np.concatenate(
            (efeature_v2f, hop.reshape(48, 6, 1)), axis=2)

        return hop, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f

    def get_high_factor_structure(self, y):
        hop = self.get_highorder_feature(y)

        nn_idx_node = self.nn_idx[:96, :3]
        # print(nn_idx_node)
        feature_h = np.take(hop, nn_idx_node.reshape(-1) - 96,
                            axis=0).reshape(96, 3, 6).astype(np.float32)

        efeature_node = np.concatenate(
            [feature_h, np.repeat(y.reshape(96, 1, 1), 3, axis=1)], axis=2)
        efeature_node_pad = np.zeros_like(efeature_node).astype(np.float32)

        efeature_node = np.concatenate(
            (efeature_node, efeature_node_pad), axis=1)

        efeature_hop = np.repeat(hop.reshape(48, 1, 6), 6, axis=1)
        efeature_hop = np.concatenate(
            (efeature_hop, hop.reshape(48, 6, 1)), axis=2)

        efeature = np.concatenate((efeature_node, efeature_hop), axis=0)

        return self.nn_idx, self.etype, efeature, hop


class Codes(Dataset):
    def __init__(self, filename, train=True):
        self.data = torch.load(filename)

        self.generator = ldpc_graph_structure_generator()

    def __len__(self):
        return len(self.data['noizy_sg'])

    def __getitem__(self, idx):

        noizy_sg = self.data['noizy_sg'][idx].numpy()
        orig = self.data['gts'][idx].numpy()

        # orig, trans, trans_noizy = gen_data_item(snr_db, sigma_b, train=False)
        hop, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f = self.generator.get_mpnn_sp_structure(
            noizy_sg)
        node_feature = torch.cat([self.data['noizy_sg'][idx].view(-1, 1),
                                  self.data['snr_dbs'][idx].view(-1, 1)], dim=1).numpy().T
        hop_feature = np.expand_dims(hop.T, -1).astype(np.float32)
        efeature_v2f = np.transpose(efeature_v2f, [2, 0, 1]).astype(np.float32)
        efeature_f2v = np.transpose(efeature_f2v, [2, 0, 1]).astype(np.float32)
        sigma_b = self.data['sigma_b'][idx].item()

        return node_feature, hop_feature, nn_idx_f2v.astype(np.int), nn_idx_v2f.astype(np.int), efeature_f2v, efeature_v2f, orig.astype(np.int), float(sigma_b)


class ContinusCodesBasic(Dataset):
    def __init__(self, lenth=10000, sigma_b=0, snr_db=0, burst_prob=0.05):
        self.len = lenth
        self.sigma_b = sigma_b
        self.snr_db = snr_db
        self.generator = ldpc_graph_structure_generator()
        self.burst_prob = burst_prob

    def __len__(self):
        return self.len

    def generate_feature(self, sigma_b, snr_db):
        sigma_b = self.sigma_b
        snr_db = self.snr_db

        orig, trans, trans_noizy, _ = gen_data_item(
            snr_db, sigma_b, burst_prob=self.burst_prob)
        nn_idx, etype, efeature, hop = self.generator.get_high_factor_structure(
            trans_noizy)
        # hop = self.generator.get_highorder_feature(trans_noizy)

        node_feature = np.asarray([[elem, snr_db]
                                   for elem in trans_noizy], dtype=np.float32).T
        hop_feature = np.expand_dims(hop.T, -1)
        etype = etype.astype(np.float32)
        efeature = np.transpose(efeature, [2, 0, 1]).astype(np.float32)
        return node_feature.astype(np.float32), hop_feature.astype(np.float32), nn_idx, etype, efeature, trans.astype(np.int), sigma_b

    def __getitem__(self, index):
        return self.generate_feature(self.sigma_b, self.snr_db)


class ContinusCodes(ContinusCodesBasic):
    def __init__(self,
                 length=10000,
                 sigma_b=[0, 1, 2, 3, 4, 5],
                 snr_db=[0, 1, 2, 3, 4],
                 burst_prob=0.05):
        super(ContinusCodes, self).__init__(
            length, sigma_b, snr_db, burst_prob)

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        sigma_b = np.random.choice(self.sigma_b)
        snr_db = np.random.choice(self.snr_db)
        return self.generate_feature(sigma_b, snr_db)


class ContinousCodesSP(Dataset):
    def __init__(self, snr=None):
        self.len = 10000
        self.sigma_b = [0, 1, 2, 3, 4, 5]
        if snr is not None:
            self.snr_db = [snr]
        else:
            self.snr_db = [0, 1, 2, 3, 4]
        self.generator = ldpc_graph_structure_generator()

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        sigma_b = np.random.choice(self.sigma_b)
        snr_db = np.random.choice(self.snr_db)

        trans_noizy, _, orig = gen_data_item(snr_db, sigma_b, train=False)
        # orig, trans, trans_noizy = gen_data_item(snr_db, sigma_b, train=False)
        hop, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f = self.generator.get_mpnn_sp_structure(
            trans_noizy)
        node_feature = np.asarray([[elem, snr_db]
                                   for elem in trans_noizy], dtype=np.float32).T
        hop_feature = np.expand_dims(hop.T, -1).astype(np.float32)
        efeature_v2f = np.transpose(efeature_v2f, [2, 0, 1]).astype(np.float32)
        efeature_f2v = np.transpose(efeature_f2v, [2, 0, 1]).astype(np.float32)

        return node_feature, hop_feature, nn_idx_f2v.astype(np.int), nn_idx_v2f.astype(np.int), efeature_f2v, efeature_v2f, orig.astype(np.int), float(sigma_b)
