import torch
from .mp_nn import mp_conv_v2
from .mp_nn_residual import mp_conv_residual
from .pooling import global_pooling
from .base_model import base_mp_nn


class mp_sequential(base_mp_nn):
    """
    Sequential module to handle a list of modules.
    It is almost the same as the torch.nn.Sequential.
    """

    def __init__(self, *module_list):
        super(mp_sequential, self).__init__()
        self.module_list = []
        for midx, mod in enumerate(module_list):
            self.add_module(str(midx), mod)
            self.module_list.append(mod)

    def forward(self, node_feature, *argv):
        """
        :param node_feature: features on n nodes, organized as b x c x n x 1
        :param argv: graph and edge features.
        """
        extra_res = []
        for m in self.module_list:
            if isinstance(m, base_mp_nn):
                node_feature = m(node_feature, *argv)
            else:
                node_feature = m(node_feature)
            if isinstance(node_feature, tuple):
                extra_feature = list(node_feature[1:])
                node_feature = node_feature[0]
                extra_res += extra_feature
        if len(extra_res) == 0:
            return node_feature
        else:
            return node_feature, extra_res
