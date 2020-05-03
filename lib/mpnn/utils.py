import torch


def pairwise_distance(node_feature):
    batch_size = node_feature.shape[0]
    node_feature = node_feature.squeeze()
    if batch_size == 1:
        node_feature = node_feature.unsqueeze(0)
    # print(node_feature.shape)
    assert (len(node_feature.shape) == 3)

    node_feature_t = node_feature.permute(0, 2, 1)
    node_feature_inner = -2 * torch.bmm(node_feature_t, node_feature)
    node_feature_square = node_feature**2
    node_feature_square = node_feature_square.sum(dim=1, keepdim=True)
    node_feature_square_t = node_feature_square.permute(0, 2, 1)
    res = node_feature_square + node_feature_square_t + node_feature_inner
    return res


def get_nn_node_feature(node_feature, nn_idx):
    batch_size = nn_idx.shape[0]
    node_feature = node_feature.squeeze()
    if batch_size == 1:
        node_feature = node_feature.unsqueeze(0)
    # print(node_feature.shape)
    # print(nn_idx.shape)
    assert (batch_size == node_feature.shape[0])
    npts = nn_idx.shape[1]
    # assert (npts == nodel_feature.shape[2])
    k = nn_idx.shape[2]
    nidx = nn_idx.view(batch_size,
                       -1).unsqueeze(1).repeat(1, node_feature.shape[1], 1)
    pts_knn = node_feature.gather(2, nidx).view(batch_size, -1, npts, k)
    return pts_knn


def get_edge_feature(pts_knn, pts):
    return pts.view(pts.shape[0], pts.shape[1], pts.shape[2], 1) - pts_knn
