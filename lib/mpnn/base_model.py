import torch


class base_mp_nn(torch.nn.Module):

    NO_EXTENSION = 0
    ORIG_WITH_NEIGHBOR = 1
    ORIG_WITH_DIFF = 2

    """
    Base class for message passing neural network
    """

    def __init__(self):
        super(base_mp_nn, self).__init__()
        self.is_mp_nn = True


class max_pool_layer(torch.nn.Module):
    """
    Max pooling over all nodes
    """

    def __init__(self, dim=2):
        super(max_pool_layer, self).__init__()

        self.dim = dim

    def forward(self, input):
        res, _ = input.max(dim=self.dim, keepdim=True)
        return res


class flatten(torch.nn.Module):
    """
    Flatten the feature for FC layers. 
    """

    def forward(self, input):
        return input.view(input.size(0), -1)


class iid_mapping(torch.nn.Module):
    """
    i.i.d mapping for each node
    """

    def __init__(self, nin, nout, bias=True):
        """
        :param nin: number of input units
        :param nout: number of output units
        :param bias: use bias or not
        """
        super(iid_mapping, self).__init__()
        self.main = torch.nn.Sequential(
            torch.nn.Conv2d(nin, nout, 1, bias=bias),
            torch.nn.ReLU(inplace=True))

    def forward(self, x):
        return self.main(x)


class iid_mapping_bn(torch.nn.Module):
    """
    i.i.d mapping with batch norm for each node
    """

    def __init__(self, nin, nout, bias=True, bn=True):
        """
        :param nin: number of input units
        :param nout: number of output units
        :param bias: use bias or not
        """
        super(iid_mapping_bn, self).__init__()
        self.main = torch.nn.Sequential(
            torch.nn.Conv2d(nin, nout, 1, bias=bias),
            torch.nn.BatchNorm2d(nout), torch.nn.ReLU(inplace=True))

    def forward(self, x):
        return self.main(x)


class iid_mapping_in(torch.nn.Module):
    def __init__(self, nin, nout, bias=True):
        super(iid_mapping_in, self).__init__()
        self.main = torch.nn.Sequential(torch.nn.Conv2d(nin, nout, 1, bias=bias),
                                        torch.nn.InstanceNorm2d(nout), torch.nn.ReLU())

    def forward(self, x):
        return self.main(x)
