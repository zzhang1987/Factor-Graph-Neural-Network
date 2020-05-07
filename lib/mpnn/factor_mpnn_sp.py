import torch
from .mp_nn import mp_conv_v2, mp_conv_type
from .base_model import iid_mapping, iid_mapping_bn, iid_mapping_in
from .mp_nn_residual import mp_conv_residual
from .base_model import base_mp_nn


class FVModuleBase (torch.nn.Module):
    def __init__(self):
        super(FVModuleBase, self).__init__()


class FVModule(FVModuleBase):
    def __init__(self, nin, nout, nedge_types, with_bn=True):
        super(FVModule, self).__init__()

        self.main_module = mp_conv_v2(
            nin, nout, nedge_types, bn=with_bn, extension=mp_conv_type.NO_EXTENSION, aggregtor='max')

    def forward(self, n_or_f_feature, nn_idx, etype):
        return self.main_module(n_or_f_feature, nn_idx, etype)


class FactorNN(torch.nn.Module):

    def __init__(self, node_feature_dim,
                 factor_feature_dim_list,
                 dim_mapping_list,
                 netype_list,
                 nclass=2,
                 gnn_immediate_dim=64,
                 max_mpnn_dim=64,
                 final_filter=None,
                 skip_link={}):
        super(FactorNN, torch.nn.Module).__init__()

        self.node_feature_dim = node_feature_dim
        self.map_dim = dim_mapping_list[0]

        self.node_mapping_module = iid_mapping(
            self.node_feature_dim, self.map_dim)

        self.nfactor_types = len(factor_feature_dim_list)

        self.factor_mapping_modules = []

        for dim in factor_feature_dim_list:
            self.factor_mapping_modules.append(iid_mapping(dim, self.map_dim))

        for idx, m in enumerate(self.factor_mapping_module):
            self.add_module('factor_mapping_modules_{}'.format(idx), m)

        self.f2v_modules = []
        self.v2f_modules = []
        self.f2f_modules = []
        self.v2v_modules = []

        for idx in range(len(dim_mapping_list) - 1):
            nin = dim_mapping_list[idx]
            nout = dim_mapping_list[idx + 1]

            self.v2v_modules.append(iid_mapping_in(nin, nout))
            self.add_module('v2v_{}'.format(idx), self.v2v_modules[-1])

            cf2v_module = []
            cv2f_module = []
            cf2f_module = []
            for jidx in range(self.nfactor_types):
                netype = netype_list[jidx]

                cf2f_module.append(iid_mapping_in(nin, nout))

                if nin == nout:
                    cf2v_module.append(mp_conv_residual(
                        nin, gnn_immediate_dim, netype, extension=mp_conv_type.NO_EXTENSION))
                    cv2f_module.append(mp_conv_residual(
                        nin, gnn_immediate_dim, netype, extension=mp_conv_type.NO_EXTENSION))

                elif nin <= max_mpnn_dim and nout <= max_mpnn_dim:
                    cf2v_module.append(mp_conv_v2(
                        nin, nout, netype, extension=mp_conv_type.NO_EXTENSION))
                    cv2f_module.append(mp_conv_v2(
                        nin, nout, netype, extension=mp_conv_type.NO_EXTENSION))
                else:
                    cf2v_module.append(iid_mapping_in(nin, nout))
                    cv2f_module.append(iid_mapping_in(nin, nout))

                self.add_module('f2v_{}_{}'.format(idx, jidx), cf2v_module[-1])
                self.add_module('v2f_{}_{}'.format(idx, jidx), cv2f_module[-1])
                self.add_module('f2f_{}_{}'.format(idx, jidx), cf2f_module[-1])

            self.f2f_modules.append(cf2f_module)
            self.f2v_modules.append(cf2v_module)
            self.v2f_modules.append(cv2f_module)
            self.skip_link = skip_link

        final_dim = nclass if nclass > 2 else 1

        self.final_classifier = torch.nn.Sequential(
            torch.nn.Conv2d(dim_mapping_list[-1], 128, 1),
            torch.nn.InstanceNorm2d(128),
            torch.nn.ReLU(inplace=True),
            torch.nn.Conv2d(128, final_dim, 1, bias=True))

    def forward(self, node_feature: torch.Tensor,
                hop_features: list,
                nn_idx_f2v: list,
                nn_idx_v2f: list,
                etype_f2v: list,
                etype_v2f: list):

        nnode_feature = self.node_mapping_module(node_feature)
        nhop_feature = [m(f) for f, m in zip(
            hop_features, self.factor_mapping_modules)]

        all_inter_features = []

        for idx in range(len(self.v2f_modules)):

            nfeature = self.v2v_modules[idx](nnode_feature)
            nffeature = [m(f)
                         for f, m in zip(nhop_feature, self.f2f_modules[idx])]
            for jidx in range(len(self.f2v_modules[idx])):

                nv = self.f2v_modules[idx][jidx](
                    nhop_feature[jidx], nn_idx_f2v[jidx], etype_f2v[jidx])
                nfeature = nfeature + nv

                nf = self.v2f_modules[idx][jidx](
                    nnode_feature, nn_idx_v2f[jidx], etype_v2f[jidx])
                nffeature[jidx] = nffeature[jidx][nf]

            all_inter_features.append([nfeature, nffeature])
            nnode_feature = nfeature
            nhop_feature = nffeature

        final_res = self.final_classifier(nhop_feature)

        return final_res
