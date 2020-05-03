import torch
from .base_model import base_mp_nn


def add_agg(*input):
    for idx, pi in enumerate(input):
        if idx == 0:
            res = pi
        else:
            res = res + pi
    return res


class parallel_net(base_mp_nn):
    """
    A module to handle several parallel neural networks. The input 
    feature will be put into n parallel neural networks to produce 
    n output features. Then an aggregator will be used to aggregate
    information from the n output features to produce a single output
    feature. 
    """

    def __init__(self, *module_list, aggregator=add_agg):
        """
        :param module_list: a list of torch modules.
        :param aggregator: a aggregator function. 
        """
        super(parallel_net, self).__init__()
        self.aggregator = aggregator
        self.module_list = []
        for midx, mod in enumerate(module_list):
            self.add_module(str(midx), mod)
            self.module_list.append(mod)

    def forward(self, node_feature, nn_idx=None, etype=None):
        res_list = []
        for m in self.module_list:
            if isinstance(m, base_mp_nn):
                res_list.append(m(node_feature, nn_idx, etype))
            else:
                res_list.append(m(node_feature))

        return self.aggregator(*res_list)
