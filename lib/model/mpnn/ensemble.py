import torch
from .mp_nn import mp_conv_v2
from .mp_nn_residual import mp_conv_residual
from .pooling import global_pooling
from .base_model import base_mp_nn


class mp_ensemble(base_mp_nn):
    def __init__(self, model1, model2, model3):
        super(mp_ensemble, self).__init__()
        self.model1 = model1
        self.model2 = model2
        self.model3 = model3

    def forward(self, node_feature, nn_idx, etype, *argv):
        x1 = self.model1(node_feature, nn_idx, etype)
        x2 = self.model2(node_feature, *argv)
        x = torch.cat((x1, x2), dim=1)
        return self.model3(x)