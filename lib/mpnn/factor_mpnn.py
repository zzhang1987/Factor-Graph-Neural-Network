import torch
from .mp_nn import mp_conv_v2
from .base_model import iid_mapping, iid_mapping_bn
from .mp_nn_residual import mp_conv_residual
from .base_model import base_mp_nn


class factor_mpnn(torch.nn.Module):
    """
    Factor GNN class. This module is used to passing messages between factors and nodes. 
    """

    def __init__(self,
                 node_feature_dim,
                 factor_feature_dim_list,
                 dim_mapping_list,
                 netype_list,
                 gnn_immediate_dim=64,
                 max_mpnn_dim=64,
                 final_filter=None,
                 skip_link={}):
        """
        :param node_feature_dim: dimension of input node feature. 
        :param factor_feature_dim_list: a list of factor feature dimensions.
        :param dim_mapping_list: a list of feature dimensions
        :param netype_list: a list of integers, a hyper parameter for Q model.
        :param gnn_immediate_dim: number of outputs for a GNN layer
        :param max_mpnn_dim: max dimension for a GNN layer (mainly for memory constraints)
        """
        super(factor_mpnn, self).__init__()
        self.node_feature_dim = node_feature_dim
        self.map_dim = dim_mapping_list[0]

        self.mapping_modules = [
            iid_mapping(self.node_feature_dim, self.map_dim)
        ]

        self.nfactor_types = len(factor_feature_dim_list)

        for dim in factor_feature_dim_list:
            self.mapping_modules.append(iid_mapping(dim, self.map_dim))

        for idx, m in enumerate(self.mapping_modules):
            self.add_module('mapping_modules_{}'.format(idx), m)

        self.mp_nn_modules = []
        self.mp_merge_modules = []
        self.final_filter = final_filter
        for idx in range(len(dim_mapping_list) - 1):
            cmodule = []
            nin = dim_mapping_list[idx]
            nout = dim_mapping_list[idx + 1]

            for jdx in range(self.nfactor_types):

                if nin == nout:
                    cmp_nn = mp_conv_residual(nin, gnn_immediate_dim,
                                              netype_list[jdx])
                else:
                    if nin <= max_mpnn_dim and nout <= max_mpnn_dim:
                        cmp_nn = mp_conv_v2(nin, nout, netype_list[jdx])
                    else:
                        cmp_nn = torch.nn.Sequential(
                            torch.nn.Conv2d(nin, nout, 1),
                            torch.nn.InstanceNorm2d(nout),
                            torch.nn.ReLU(inplace=True))

                self.add_module('mp_nn_{}_{}'.format(idx, jdx), cmp_nn)
                cmodule.append(cmp_nn)
            self.mp_nn_modules.append(cmodule)
            if(idx < len(dim_mapping_list) - 2):
                merge_module = iid_mapping_bn(nout * self.nfactor_types, nout)
            else:
                merge_module = torch.nn.Sequential(
                    torch.nn.Conv2d(
                        nout * self.nfactor_types, 256, 1, bias=True),
                    torch.nn.BatchNorm2d(256),
                    torch.nn.LeakyReLU(),
                    torch.nn.Conv2d(256, 256, 1, bias=True),
                    torch.nn.LeakyReLU(),
                    torch.nn.Conv2d(256, nout, 1, bias=True)
                )

            self.add_module('merge_module_{}'.format(idx), merge_module)
            self.mp_merge_modules.append(merge_module)
            self.skip_link = skip_link

    def forward(self, node_features,  factor_features, graph_structures):
        """
        :param node_features: features on nodes 
        :param factor_features: list of features on different types of factors 
        :param graph_structures: factor graph structure.
        """
        nnode = node_features.shape[2]
        nfeatures = self.mapping_modules[0](node_features)

        ffeatures = [
            m(f) for f, m in zip(factor_features, self.mapping_modules[1:])
        ]

        all_inter_features = []
        for midx, modules in enumerate(self.mp_nn_modules):
            cnfeatures = []
            cffeatures = []
            for fidx, (f, m) in enumerate(zip(ffeatures, modules)):
                cfeature = torch.cat([nfeatures, f], dim=2)
                nn_idx, etype = graph_structures[fidx]
                if isinstance(m, base_mp_nn):
                    cfeature = m(cfeature.contiguous(), nn_idx, etype)
                else:
                    cfeature = m(cfeature.contiguous())

                cnfeatures.append(cfeature[:, :, :nnode, :])
                cffeatures.append(cfeature[:, :, nnode:, :])
            cnfeatures = torch.cat(cnfeatures, dim=1)
            nfeatures = self.mp_merge_modules[midx](cnfeatures)
            # print('Layer _{}'.format(midx), self.mp_merge_modules[midx])
            ffeatures = cffeatures

            if midx in self.skip_link.keys():
                tidx = self.skip_link[midx]
                onfeatures, offeatures = all_inter_features[tidx]
                nfeatures = nfeatures + onfeatures
                nffeatures = [
                    ff + off for (ff, off) in zip(ffeatures, offeatures)]
                ffeatures = nffeatures
            all_inter_features.append([nfeatures, ffeatures])

        if self.final_filter is not None:
            nfeatures = self.final_filter(nfeatures, node_features)
        # print(nfeatures[0, ...].squeeze())
        # print(nfeatures[1, ...].squeeze())
        return nfeatures, ffeatures
