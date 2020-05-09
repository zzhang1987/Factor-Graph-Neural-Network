import torch
from torch.utils.data import Dataset

import numpy as np
import pathlib
import os
import tempfile
import subprocess
from .MNC import init_seed, s2t, t2y
import pdb


class LDPCGenerator:
    def __init__(self):
        self.code_len = 48
        self.Gfile = os.path.join(
            pathlib.Path(__file__).parent.absolute(),
            './MNC/codes/96.3.963/G'
        )
        self.ldpc_file = os.path.join(
            pathlib.Path(__file__).parent.absolute(),
            '../../ldpc_data/96.3.963'
        )
        self.get_factor_sturcture()

    def get_factor_sturcture(self):
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

    def __call__(self, snr_db, sigma_b):

        s = np.random.randint(0, 2, self.code_len)
        t = s2t(s, 48, 48, self.Gfile, True)
        y = t2y(t, snr_db, sigma_b, 0.05)
        return s, t, y

    def get_highorder_feature(self, y):

        return np.take(y, self.factors.reshape(-1)).reshape(self.factors.shape)

    def get_mpnn_sp_structure(self, y):
        hop = self.get_highorder_feature(y)
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
        self.node_feature = self.data['node_feature']
        self.hop_feature = self.data['hop_feature']
        self.y = self.data['y']

        if not train:
            self.sigma_b = self.data['sigma_b']
        else:
            self.sigma_b = torch.zeros((len(self.node_feature)))

        self.len = len(self.node_feature)
        self.generator = LDPCGenerator()

    def __len__(self):
        return self.len

    def __getitem__(self, idx):

        node_feature = self.node_feature[idx].numpy().squeeze().T[:, 0]
        nn_idx, etype, efeature, hop = self.generator.get_high_factor_structure(
            node_feature)
        # hop = self.generator.get_highorder_feature(node_feature)
        hop_feature = np.expand_dims(hop.T, -1)
        etype = etype.astype(np.float32)
        efeature = np.transpose(efeature, [2, 0, 1]).astype(np.float32)
        node_feature = self.node_feature[idx].squeeze()
        # print(node_feature[1, :])
        # node_feature[1, :] = 10 * torch.log10(node_feature[1, :])
        # print(node_feature)
        # pdb.set_trace()
        return node_feature.unsqueeze(-1), hop_feature, nn_idx, etype, efeature, self.y[idx], self.sigma_b[idx]


class Codes_SP(Codes):
    def __init__(self, filename, train=True):
        super(Codes_SP, self).__init__(filename, train)

    def __getitem__(self, idx):
        node_feature = self.node_feature[idx].numpy().squeeze().T[:, 0]
        hop, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f = self.generator.get_mpnn_sp_structure(
            node_feature)
        node_feature = self.node_feature[idx].squeeze()

        hop_feature = np.expand_dims(hop.T, -1).astype(np.float32)
        efeature_v2f = np.transpose(efeature_v2f, [2, 0, 1]).astype(np.float32)
        efeature_f2v = np.transpose(efeature_f2v, [2, 0, 1]).astype(np.float32)

        return node_feature, hop_feature, nn_idx_f2v.astype(np.int), nn_idx_v2f.astype(np.int), efeature_f2v, efeature_v2f, self.y[idx], self.sigma_b[idx]


class ContinusCodes(Dataset):
    def __init__(self):
        self.len = 10000
        self.sigma_b = [0, 1, 2, 3, 4, 5]
        self.snr_db = [0, 1, 2, 3, 4]
        self.generator = LDPCGenerator()

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        sigma_b = np.random.choice(self.sigma_b)
        snr_db = np.random.choice(self.snr_db)

        orig, trans, trans_noizy = self.generator(snr_db, sigma_b)
        nn_idx, etype, efeature, hop = self.generator.get_high_factor_structure(
            trans_noizy)
        # hop = self.generator.get_highorder_feature(trans_noizy)

        node_feature = np.asarray([[elem, snr_db]
                                   for elem in trans_noizy], dtype=np.float32).T
        hop_feature = np.expand_dims(hop.T, -1)
        etype = etype.astype(np.float32)
        efeature = np.transpose(efeature, [2, 0, 1]).astype(np.float32)
        return node_feature.astype(np.float32), hop_feature.astype(np.float32), nn_idx, etype, efeature, trans.astype(np.int), sigma_b


class ContinousCodesSP(Dataset):
    def __init__(self):
        self.len = 10000
        self.sigma_b = [0, 1, 2, 3, 4, 5]
        self.snr_db = [0, 1, 2, 3, 4]
        self.generator = LDPCGenerator()

    def __len__(self):
        return self.len

    def __getitem__(self, index):
        sigma_b = np.random.choice(self.sigma_b)
        snr_db = np.random.choice(self.snr_db)

        orig, trans, trans_noizy = self.generator(snr_db, sigma_b)
        hop, nn_idx_f2v, nn_idx_v2f, efeature_f2v, efeature_v2f = self.generator.get_mpnn_sp_structure(
            trans_noizy)
        node_feature = np.asarray([[elem, snr_db]
                                   for elem in trans_noizy], dtype=np.float32).T
        hop_feature = np.expand_dims(hop.T, -1).astype(np.float32)
        efeature_v2f = np.transpose(efeature_v2f, [2, 0, 1]).astype(np.float32)
        efeature_f2v = np.transpose(efeature_f2v, [2, 0, 1]).astype(np.float32)

        return node_feature, hop_feature, nn_idx_f2v.astype(np.int), nn_idx_v2f.astype(np.int), efeature_f2v, efeature_v2f, trans.astype(np.int), float(sigma_b)
