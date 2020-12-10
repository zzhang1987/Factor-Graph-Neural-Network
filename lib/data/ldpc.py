import numpy as np
from .MNC import s2t, t2y, y2b, zb2x, init_seed
import os
import logging



def gen_data_item(snr_db, sigma_b, burst_prob=0.05, train=True):
    
    s = np.random.randint(0, 2, 48)

    # Gfile for coding
    Gfile = os.path.join(os.path.dirname(__file__),
                         '../../ldpc_codes/96.3.963/G')
    # encode
    t = s2t(s, 48, 48, Gfile, True)
    # add noise
    y = t2y(t, snr_db, sigma_b, burst_prob)
    if train:
        z = y2b(y, snr_db)
        Afile = os.path.join(os.path.dirname(__file__),
                             '../../ldpc_codes/96.3.963/A2')
        # decode using sum-product
        x = zb2x(z, 48, 48, Afile, 1, 100)
        error = np.sum(np.reshape(x, [48]) != np.reshape(s, [48])) / 48
        return y, np.ones((9, 48)), t, error
    else:
        return y, np.ones((9, 48)), t
