import torch
import logging
from .base_model import base_mp_nn


def max_pool(feature, dim):
    res, _ = feature.max(dim=dim, keepdim=True)
    return res


class global_pooling(base_mp_nn):
    """
    The global pooling layer.
    """

    def __init__(self,
                 orig_mapper=None,
                 gfeature_mapper=None,
                 pool_func=lambda x: max_pool(x, dim=2)):
        super(global_pooling, self).__init__()
        self.orig_mapper = orig_mapper
        self.gfeature_mapper = gfeature_mapper
        self.pool_func = pool_func
        if orig_mapper is not None:
            self.add_module('orig_mapper', self.orig_mapper)
        if gfeature_mapper is not None:
            self.add_module('gfeature_mapper', self.gfeature_mapper)

    def forward(self, node_feature, nn_idx, etype):
        nnodes = node_feature.shape[2]

        logging.debug('node_feature shape {}'.format(node_feature.shape))

        gfeature = self.pool_func(node_feature)
        logging.debug('gfeature shape {}'.format(gfeature.shape))

        if self.orig_mapper is not None:
            node_feature = self.orig_mapper(node_feature, nn_idx, etype)
        if self.gfeature_mapper is not None:
            gfeature = self.gfeature_mapper(gfeature).repeat(1, 1, nnodes, 1)
        else:
            gfeature = gfeature.repeat(1, 1, nnodes, 1)
        logging.debug('gfeature shape {}'.format(gfeature.shape))
        node_feature = torch.cat([node_feature, gfeature], dim=1)
        return node_feature
