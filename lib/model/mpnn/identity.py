import torch


class identity_module(torch.nn.Module):
    """
    Identity module, returns the input
    """

    def __init__(self):
        super(identity_module, self).__init__()

    def forward(self, input):
        return input
