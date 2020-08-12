from torch.utils.data import Dataset
try:
    import cPickle as pickle
except:
    import pickle
import time
import numpy as np
from ad3 import factor_graph as fg
import os


class RandomPGMHopVarChainsze(Dataset):
    def __init__(self, chain_length_range, hop_order=9, ret_efeature_pw=True, size=10000, save_dir='./train_syn_var_size'):
        self.chain_length_range = chain_length_range
        self.hop_order = hop_order if hop_order >> 1 else hop_order + 1
        self.half_hop = self.hop_order >> 1
        self.ret_efeature_pw = ret_efeature_pw
        self.size = size

        if not os.path.exists(save_dir):
            os.mkdir(save_dir)

        self.save_dir = save_dir

    def __len__(self):
        return self.size

    def _generate_graph(self, chain_length, lops, pws, cap):
        g = fg.PFactorGraph()
        var_list = []
        for i in range(chain_length):
            v = g.create_multi_variable(2)
            v.set_log_potentials(lops[:, i].astype(np.float))
            var_list.append(v)

        for i in range(chain_length - 1):
            g.create_factor_dense(
                [var_list[i], var_list[i + 1]], pws[i][0].astype(np.float))

        for i in range(chain_length - self.hop_order + 1):
            v_list = [
                var_list[j].get_state(1) for j in range(i, i + self.hop_order)
            ]
            g.create_factor_budget(v_list, cap[i + self.half_hop])
        return g

    def _get_node_feature(self, chain_length):
        lops = np.random.uniform(0.0, 1.0, (chain_length, 2))
        return np.transpose(lops.astype(np.float32), [1, 0])

    def _get_edge_feature_pw(self, chain_length):
        pws = np.zeros(shape=[chain_length, 2, 4], dtype=np.float32)
        for i in range(chain_length - 1):
            # pws_to_right = np.random.randn(2, 2)
            pws_to_right = np.zeros([2, 2])
            pws_to_right[1, 1] = np.random.uniform(0, 2)
            pws_to_left = np.transpose(pws_to_right)
            pws[i] = [list(pws_to_right.reshape(-1)),
                      list(pws_to_left.reshape(-1))]

        efeature = np.zeros(shape=[chain_length, 3, 4], dtype=np.float32)
        for i in range(chain_length):
            e_self = np.zeros(4)
            e_left = pws[i-1][1] if i > 0 else e_self.copy()
            e_right = pws[i][0] if i < chain_length - \
                1 else e_self.copy()
            efeature[i, 0] = e_left
            efeature[i, 1] = e_self
            efeature[i, 2] = e_right
        return pws, np.transpose(efeature, [2, 0, 1])

    def _generate_edge_feature_hop(self, chain_length):
        cap = list(np.random.randint(
            low=1, high=self.hop_order, size=chain_length))
        half_hop = self.hop_order >> 1

        max_cap = np.zeros(self.hop_order)
        max_cap[self.hop_order-1] = 1

        efeature = np.zeros(
            shape=[chain_length, self.hop_order], dtype=np.float32)
        for i in range(half_hop, chain_length - half_hop):
            efeature[i, cap[i]] = 1
        for i in range(half_hop):
            efeature[i, self.hop_order-1] = 1
        for i in range(chain_length - half_hop, chain_length):
            efeature[i, self.hop_order-1] = 1

        return cap, np.expand_dims(np.transpose(efeature, [1, 0]), -1)

    def __getitem__(self, index):
        save_path = os.path.join(self.save_dir, '{}.data'.format(index))
        try:
            if os.path.exists(save_path):
                with open(save_path, "rb") as f:
                    cdata = pickle.load(f)
                    return cdata
        except:
            pass

        cl_start, cl_end = self.chain_length_range
        chain_length = np.random.choice(list(range(cl_start, cl_end)))
        node_feature = self._get_node_feature(chain_length)
        pws, efeature_pw = self._get_edge_feature_pw(chain_length)
        cap, efeature_hop = self._generate_edge_feature_hop(chain_length)

        ''' exact solution '''
        g = self._generate_graph(chain_length, node_feature, pws, cap)

        val, post, _, stat = g.solve(tol=1e-6, branch_and_bound=True)
        post = np.reshape(np.asarray(post), [chain_length, 2])
        assign = np.argmax(post, axis=1)

        ''' approx solution '''
        g = self._generate_graph(chain_length, node_feature, pws, cap)

        val1, post1, _, status = g.solve(branch_and_bound=False)
        post1 = np.reshape(np.asarray(post1), [chain_length, 2])
        assign1 = np.argmax(post1, axis=1)

        if self.ret_efeature_pw:
            with open(save_path, "wb") as f:
                pickle.dump((node_feature, efeature_pw,
                             efeature_hop, assign, assign1), f)
            return node_feature, efeature_pw, efeature_hop, assign, assign1
        else:
            pws = np.expand_dims(np.transpose(pws[:, 0, :], [1, 0]), -1)
            with open(save_path, "wb") as f:
                pickle.dump(
                    (node_feature, pws, efeature_hop, assign, assign1), f)
            return node_feature, pws, efeature_hop, assign, assign1


class RandomPGMData(Dataset):
    def __init__(self, filename, pgm_type, size):
        f = open(filename, "rb")
        self.data = []
        for _ in range(size):
            self.data.append(pickle.load(f))
        self.pgm_type = pgm_type
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        if self.pgm_type == "raw":
            node_feature, assign, assign1 = self.data[idx]
            return node_feature, assign, assign1
        elif self.pgm_type == "pws":
            node_feature, edge_feature, assign, assign1 = self.data[idx]
            return node_feature, edge_feature, assign, assign1
        elif self.pgm_type == "hops":
            node_feature, efeature_pw, efeature_hop, assign, assign1 = self.data[idx]
            # print("res", node_feature, efeature_pw,
            #       efeature_hop, assign, assign1)
            return node_feature, efeature_pw, efeature_hop, assign, assign1
        else:
            print("pgm type error")
            exit(-1)


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    np.random.seed(((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) +
                   ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24))
