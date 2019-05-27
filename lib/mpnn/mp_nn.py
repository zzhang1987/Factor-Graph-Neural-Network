import torch
from enum import Enum
from .base_model import base_mp_nn
SyncBatchNorm = torch.nn.BatchNorm2d


class mp_conv_type(Enum):
    NO_EXTENSION = 0
    ORIG_WITH_NEIGHBOR = 1
    ORIG_WITH_DIFF = 2


class mp_conv_v2(base_mp_nn):
    """
    Message passing neural network.
    """

    def __init__(self,
                 nin,
                 nou,
                 nedge_types,
                 bias=True,
                 bn=True,
                 extension=mp_conv_type.ORIG_WITH_DIFF,
                 activation_fn='relu',
                 aggregtor='max'):
        """
        :param nin: dimension of input features
        :param nou: dimension of output features
        :param nedge_types: number of hidden layers in the Q net.
        :param bias: use bias or not
        :param bn: use batch norm or not
        """
        super(mp_conv_v2, self).__init__()

        self.nin = nin
        self.nou = nou
        self.nedge_types = nedge_types
        self.extension = extension

        if self.extension == mp_conv_type.NO_EXTENSION:
            self.filters = torch.nn.Parameter(
                torch.zeros(nin, nou * nedge_types, dtype=torch.float32))

        elif self.extension == mp_conv_type.ORIG_WITH_DIFF or self.extension == mp_conv_type.ORIG_WITH_NEIGHBOR:
            self.filters = torch.nn.Parameter(
                torch.zeros(2 * nin, nou * nedge_types, dtype=torch.float32))
        else:
            raise ValueError("extension must one of mp_conv_type")
        self.filters.data.uniform_(-0.01, 0.01)

        if bias:
            self.bias = torch.nn.Parameter(torch.zeros(nou))
            self.bias.data.uniform_(0, 0.05)
        else:
            self.bias = None

        if bn:
            self.bn = SyncBatchNorm(nou)
        else:
            self.bn = None
        if isinstance(activation_fn, torch.nn.Module):
            self.activation_fn = activation_fn
        elif activation_fn == 'relu':
            self.activation_fn = torch.nn.ReLU(inplace=True)
        else:
            self.activation_fn = None

        if isinstance(aggregtor, str):
            if aggregtor == 'max':

                def agg_max(x):
                    res, *_ = torch.max(x, dim=3, keepdim=True)
                    return res

                self.aggregtor = agg_max
            elif aggregtor == 'mean':
                self.aggregtor = lambda x: torch.mean(x, dim=3, keepdim=True)

        else:
            self.aggregtor = aggregtor

    def to_edge_feature(self, node_feature, nn_idx):
        batch_size = nn_idx.shape[0]
        node_feature = node_feature.squeeze()  # shape n x b x c
        if batch_size == 1:
            node_feature = node_feature.unsqueeze(0)
        # print(node_feature.shape)
        # print(nn_idx.shape)
        assert (batch_size == node_feature.shape[0])
        npts = nn_idx.shape[1]
        # assert (npts == node_feature.shape[1])
        k = nn_idx.shape[2]

        nidx = nn_idx.view(batch_size, -1).unsqueeze(2).repeat(
            1, 1, node_feature.shape[2])  # shape n x b x k

        # print(node_feature.shape)
        # print(nidx.shape)

        pts_knn = node_feature.gather(1, nidx).view(batch_size, npts, k, -1)

        return pts_knn

    def forward(self, x, nn_idx, etype):
        # print(x.shape)
        # print(self.nin, self.nou, self.nedge_types)
        batch_size = x.shape[0]
        nin = x.shape[1]
        nnodes = x.shape[2]
        k = nn_idx.shape[2]
        nedge_type = etype.permute(0, 2, 3, 1).contiguous().view(
            -1, self.nedge_types, 1)
        if self.extension == mp_conv_type.NO_EXTENSION:
            node_feature = x.permute(0, 2, 3, 1).contiguous().view(
                batch_size * nnodes, nin)
            node_feature = node_feature.matmul(
                self.filters.view(nin, self.nou * self.nedge_types)).view(
                    batch_size, nnodes, self.nou * self.nedge_types)
            edge_feature = self.to_edge_feature(node_feature, nn_idx).view(
                -1, self.nou, self.nedge_types)

            edge_feature = edge_feature.bmm(nedge_type).view(
                batch_size, nn_idx.shape[1], k, self.nou)

        else:

            node_feature = x.permute(0, 2, 3, 1).contiguous()
            # print(nn_idx.shape)
            # print(etype.shape)
            efeature = self.to_edge_feature(node_feature, nn_idx)
            # print(efeature.shape)
            if self.extension == mp_conv_type.ORIG_WITH_DIFF:
                efeature = node_feature - efeature

            node_feature = node_feature.repeat(1, 1, efeature.shape[-2], 1)
            # print(node_feature.shape)
            # print(efeature.shape)
            efeature = torch.cat([node_feature, efeature], dim=3)
            # print(node_feature.shape, self.filters[0].view(
            #     nin, self.nou * self.nedge_types).shape)
            edge_feature = efeature.view(-1, 2 * nin).matmul(
                self.filters).view(batch_size, nnodes, efeature.shape[-2],
                                   self.nou * self.nedge_types)

            edge_feature = edge_feature.view(
                -1, self.nou,
                self.nedge_types).bmm(nedge_type).view(batch_size, nnodes, k,
                                                       self.nou)
        nfeature = edge_feature.permute(0, 3, 1, 2).contiguous()

        if self.aggregtor is not None:
            nfeature = self.aggregtor(nfeature)
            # print(nfeature.shape)
        if self.bias is not None:
            # print(nfeature.shape)
            # print(self.bias.shape)
            nfeature = nfeature + self.bias.view(1, self.nou, 1, 1)
        if self.bn is not None:
            nfeature = self.bn(nfeature)

        if self.activation_fn is not None:
            nfeature = self.activation_fn(nfeature)

        return nfeature


class gconv_residual(torch.nn.Module):
    def __init__(self, nin, nmed, netype, with_residual=True, version='v2'):
        super(gconv_residual, self).__init__()
        self.conv1 = torch.nn.Sequential(torch.nn.Conv2d(nin, nmed, 1),
                                         SyncBatchNorm(nmed),
                                         torch.nn.ReLU(inplace=True))
        self.mp_conv = mp_conv_v2(nmed, nmed, netype)

        self.conv2 = torch.nn.Sequential(torch.nn.Conv2d(nmed, nin, 1),
                                         SyncBatchNorm(nin),
                                         torch.nn.ReLU(inplace=True))
        self.with_residual = with_residual

    def forward(self, node_feature, nn_idx, etype):
        nfeature = self.conv1(node_feature)
        nfeature = self.mp_conv(nfeature, nn_idx, etype)
        nfeature = self.conv2(nfeature)

        if self.with_residual:
            nfeature = nfeature + node_feature

        return nfeature
