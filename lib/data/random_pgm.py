import numpy as np
import torch
from torch.utils.data import Dataset
import numpy as np
from ad3 import factor_graph as fg

len = 100000

class RandomPGM(Dataset):
    def __init__(self, chain_length, cap, transition, hop_order=9, size=len):
        self.chain_length = chain_length
        self.cap = cap
        self.hop_order = hop_order
        self.transition = transition
        self.size = size

    def __len__(self):
        return self.size

    def __getitem__(self, index):

        lops = np.random.uniform(0.0, 1.0, (self.chain_length, 2))
        pws = self.transition

        hk = self.hop_order // 2

        g = fg.PFactorGraph()
        var_list = []
        for i in range(self.chain_length):
            v = g.create_multi_variable(2)
            v.set_log_potentials(lops[i])
            var_list.append(v)

        for i in range(self.chain_length - 1):
            g.create_factor_dense([var_list[i], var_list[i + 1]], pws)

        for i in range(self.chain_length - self.hop_order + 1):
            v_list = [
                var_list[j].get_state(1) for j in range(i, i + self.hop_order)
            ]
            g.create_factor_budget(v_list, self.cap)

        val, post, _, stat = g.solve(tol=1e-6, branch_and_bound=True)

        post = np.reshape(np.asarray(post), [self.chain_length, 2])
        assign = np.argmax(post, axis=1)

        node_feature = np.transpose(lops.astype(np.float32), [1, 0])

        g = fg.PFactorGraph()
        var_list = []
        for i in range(self.chain_length):
            v = g.create_multi_variable(2)
            v.set_log_potentials(lops[i])
            var_list.append(v)

        for i in range(self.chain_length - 1):
            g.create_factor_dense([var_list[i], var_list[i + 1]], pws)

        for i in range(self.chain_length - self.hop_order + 1):
            v_list = [
                var_list[j].get_state(1) for j in range(i, i + self.hop_order)
            ]
            g.create_factor_budget(v_list, self.cap)

        val1, post1, _, status = g.solve(branch_and_bound=False)
        post1 = np.reshape(np.asarray(post1), [self.chain_length, 2])
        assign1 = np.argmax(post1, axis=1)

        return node_feature, assign, assign1


if __name__ == '__main__':
    rpgm = RandomPGM(100, [6, 8, 3, 2, 1, 4, 2])

    assign, node_feature = rpgm[0]

    print(assign)
    print(node_feature.shape)
