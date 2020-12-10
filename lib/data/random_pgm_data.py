from torch.utils.data import Dataset
try:
    import cPickle as pickle
except:
    import pickle
import time
import numpy as np

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
            return node_feature, efeature_pw, efeature_hop, assign, assign1
        else:
            print("pgm type error")
            exit(-1)


def worker_init_fn(idx):
    t = int(time.time() * 1000.0) + idx
    np.random.seed(((t & 0xff000000) >> 24) + ((t & 0x00ff0000) >> 8) +
                ((t & 0x0000ff00) << 8) + ((t & 0x000000ff) << 24))
