import argparse
import numpy as np
from tqdm import tqdm
import torch
try:
    import cPickle as pickle
except:
    import pickle
import lib
import logging
from utils.types import str2bool


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--name',
                        type=str,
                        default="",
                        help="generated filename")
    parser.add_argument('--type', type=str, default="", help="pgm type")

    parser.add_argument('--hop',
                        type=str2bool,
                        default=False,
                        help="whether contains hop factors or not")
    parser.add_argument('--size', type=int,
                        default=100000,
                        help="The whole dataset size")

    return parser.parse_args()


def write_dataset(filename, dloader):
    with open(filename, "wb") as f:
        for dataitem in tqdm(dloader):
            batch_size = dataitem[0].shape[0]
            for i in range(batch_size):
                ndataitem = tuple([d[i].cpu().numpy() for d in dataitem])
                pickle.dump(ndataitem, f)


def generate_dataset(filename, pgm_type, hop, size):
    if pgm_type == "raw":
        # transition = list(np.random.randn(2 * 2))
        transition = [0, 0.1, 0.2, 1]
        if hop:
            rpgm = lib.data.RandomPGM(chain_length=30,
                                      cap=5,
                                      transition=transition,
                                      size=size)
            logging.info("Create Model with fixed HOP and Pairwise")
        else:
            rpgm = lib.data.RandomPGMNoHop(chain_length=30,
                                           cap=5,
                                           transition=transition,
                                           size=size)
    elif pgm_type == "pws":
        if hop:
            rpgm = lib.data.RandomPGMPw(chain_length=30,
                                        cap=5,
                                        ret_efeature=False,
                                        size=size)
        else:
            rpgm = lib.data.RandomPGMPwNoHop(chain_length=30,
                                             cap=5,
                                             ret_efeature=False,
                                             size=size)
    elif pgm_type == "hops":
        rpgm = lib.data.RandomPGMHop(chain_length=30,
                                     ret_efeature_pw=False,
                                     size=size)
    else:
        print("pgm type error")
        exit(-1)
    batch_size = 32
    dloader = torch.utils.data.DataLoader(rpgm,
                                          batch_size=batch_size,
                                          shuffle=True,
                                          num_workers=24,
                                          worker_init_fn=lib.data.worker_init_fn)

    write_dataset(filename, dloader)


if __name__ == "__main__":
    args = parse_args()
    generate_dataset(args.name, args.type, args.hop, args.size)
