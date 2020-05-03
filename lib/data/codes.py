import torch
from torch.utils.data import Dataset


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

    def __len__(self):
        return self.len

    def __getitem__(self, idx):
        return self.node_feature[idx], self.hop_feature[idx], self.y[idx], self.sigma_b[idx]