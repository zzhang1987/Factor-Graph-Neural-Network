import numpy as np
import torch
from torch.utils.data import Dataset
import numpy as np
from ad3 import factor_graph as fg

len = 100000

class RandomPGMNoHop(Dataset):
    def __init__(self, chain_length, cap, transition, hop_order=9, size=len):
        self.chain_length = chain_length
        self.cap = cap
        self.hop_order = hop_order
        self.transition = transition
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
            g.create_factor_dense([var_list[i], var_list[i + 1]], self.pws)

        return g

    def __getitem__(self, index):

        self.lops = np.random.uniform(0.0, 1.0, (self.chain_length, 2))
        self.pws = self.transition

        hk = self.hop_order // 2
        node_feature = np.transpose(self.lops.astype(np.float32), [1, 0])

        g = self._generate_graph()

        val, post, _, stat = g.solve(tol=1e-6, branch_and_bound=True)

        post = np.reshape(np.asarray(post), [self.chain_length, 2])
        assign = np.argmax(post, axis=1)


        g = self._generate_graph()
        val1, post1, _, status = g.solve(branch_and_bound=False)
        post1 = np.reshape(np.asarray(post1), [self.chain_length, 2])
        assign1 = np.argmax(post1, axis=1)

        return node_feature, assign, assign1


if __name__ == '__main__':
    rpgm = RandomPGMNoHop(100, [6, 8, 3, 2, 1, 4, 2])

    assign, node_feature = rpgm[0]

    print(assign)
    print(node_feature.shape)
