import numpy as np
import torch
from torch.utils.data import Dataset
import numpy as np
from ad3 import factor_graph as fg
try:
    import cPickle as pickle
except:
    import pickle
from tqdm import tqdm
import time

from .random_pgm_data import RandomPGMData, worker_init_fn

len = 100000

class RandomPGMHop(Dataset):
    def __init__(self, chain_length, hop_order=9, ret_efeature_pw=True, size=len):
        self.chain_length = chain_length
        self.hop_order = hop_order if hop_order >> 1 else hop_order + 1
        self.half_hop = self.hop_order >> 1
        self.ret_efeature_pw = ret_efeature_pw
        self.size = size

    def __len__(self):
        return self.size

    def _generate_graph(self):
        g = fg.PFactorGraph()
        var_list = []
        for i in range(self.chain_length):
            v = g.create_multi_variable(2)
            v.set_log_potentials(self.lops[i])
            var_list.append(v)

        for i in range(self.chain_length - 1):
            g.create_factor_dense([var_list[i], var_list[i + 1]], self.pws[i][0])

        for i in range(self.chain_length - self.hop_order + 1):
            v_list = [
                var_list[j].get_state(1) for j in range(i, i + self.hop_order)
            ]
            g.create_factor_budget(v_list, self.cap[i + self.half_hop])

        return g

    def _get_node_feature(self):
        self.lops = np.random.uniform(0.0, 1.0, (self.chain_length, 2))
        return np.transpose(self.lops.astype(np.float32), [1, 0])

    def _get_edge_feature_pw(self):
        self.pws = np.zeros(shape=[self.chain_length, 2, 4], dtype=np.float32)
        for i in range(self.chain_length - 1):
            # pws_to_right = np.random.randn(2, 2)
            pws_to_right = np.zeros([2, 2])
            pws_to_right[1, 1] = np.random.uniform(0, 2)
            pws_to_left = np.transpose(pws_to_right)
            self.pws[i] = [list(pws_to_right.reshape(-1)), list(pws_to_left.reshape(-1))]

        efeature = np.zeros(shape=[self.chain_length, 3, 4], dtype=np.float32)
        for i in range(self.chain_length):
            e_self = np.zeros(4)
            e_left = self.pws[i-1][1] if i > 0 else e_self.copy()
            e_right = self.pws[i][0] if i < self.chain_length-1 else e_self.copy()
            efeature[i, 0] = e_left
            efeature[i, 1] = e_self
            efeature[i, 2] = e_right
        return np.transpose(efeature, [2, 0, 1])

    def _generate_edge_feature_hop(self):
        self.cap = list(np.random.randint(low=1, high=self.hop_order, size=self.chain_length))
        half_hop = self.hop_order >> 1

        max_cap = np.zeros(self.hop_order)
        max_cap[self.hop_order-1] = 1

        efeature = np.zeros(shape=[self.chain_length, self.hop_order], dtype=np.float32)
        for i in range(half_hop, self.chain_length - half_hop):
            efeature[i, self.cap[i]] = 1
        for i in range(half_hop):
            efeature[i, self.hop_order-1] = 1
        for i in range(self.chain_length - half_hop, self.chain_length):
            efeature[i, self.hop_order-1] = 1

        return np.expand_dims(np.transpose(efeature, [1, 0]), -1)

        ''' passing info from node to factor '''
        efeature = np.zeros(shape=[self.chain_length, self.hop_order, self.hop_order], dtype=np.float32)
        for i in range(self.chain_length):
            cur_cap = np.zeros(self.hop_order)
            cur_cap[self.cap[i]] = 1

            for j in range(i - half_hop, i + half_hop + 1):
                efeature[i, j-i+half_hop] = max_cap.copy() \
                    if j < 0 or j >= self.chain_length else cur_cap.copy()
        return np.transpose(efeature, [2, 0, 1])

        ''' passing info from factor to node '''
        efeature2 = np.zeros(shape=[self.chain_length, self.hop_order, self.hop_order], dtype=np.float32)
        for i in range(self.chain_length):
            for j in range(i-half_hop, i+half_hop+1):
                if j < 0 or j >= self.chain_length:
                    efeature2[i, j-i+half_hop] = max_cap.copy()
                else:
                    cur_cap = np.zeros(self.hop_order)
                    cur_cap[self.cap[j]] = 1
                    efeature2[i, j-i+half_hop] = cur_cap

        return np.transpose(efeature, [2, 0, 1]), np.transpose(efeature2, [2, 0, 1])

    def __getitem__(self, index):
        node_feature = self._get_node_feature()
        efeature_pw = self._get_edge_feature_pw()
        efeature_hop = self._generate_edge_feature_hop()

        ''' exact solution '''
        g = self._generate_graph()

        val, post, _, stat = g.solve(tol=1e-6, branch_and_bound=True)

        post = np.reshape(np.asarray(post), [self.chain_length, 2])
        assign = np.argmax(post, axis=1)

        ''' approx solution '''
        g = self._generate_graph()

        val1, post1, _, status = g.solve(branch_and_bound=False)
        post1 = np.reshape(np.asarray(post1), [self.chain_length, 2])
        assign1 = np.argmax(post1, axis=1)

        if self.ret_efeature_pw:
            return node_feature, efeature_pw, efeature_hop, assign, assign1
        else:
            pws = np.expand_dims(np.transpose(self.pws[:, 0, :], [1, 0]), -1)
            return node_feature, pws, efeature_hop, assign, assign1


if __name__ == '__main__':
    rpgm = RandomPGMHop(chain_length=6, hop_order=5, ret_efeature_pw=False)

    node_feature, efeature_pw, efeature_hop, assign, assign1 = rpgm[0]

    print('node_feature', node_feature.shape, np.transpose(node_feature, [1, 0]))
    print('efeature_pw', efeature_pw.shape, np.transpose(efeature_pw, [1, 2, 0]))
    print('efeature_hop', efeature_hop.shape, np.transpose(efeature_hop, [1, 2, 0]))
    print('assign', assign)
    print('assign1', assign1)
