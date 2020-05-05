import numpy as np
import pathlib
import os
import tempfile


def get_snr(snr_db):
    return 10**(snr_db/10)


class LDPCGenerator:
    def __init__(self):
        self.code_len = 48
        self.directory = pathlib.Path(__file__).parent.absolute()

    def __call__(self, gcx, sigma):
        _, spath = tempfile.mkstemp()
        x = self.generate_orig_code(spath)
        t, y = self.transmit_cde(spath, gcx, sigma)
        os.remove(spath)
        return x, t, y

    def generate_orig_code(self, tmppath):
        s = np.random.randint(0, 2, self.code_len)
        print(tmppath)
        np.savetxt(tmppath, s, '%d')
        return s

    def decode(self, y, snr_db, x):
        gcx = get_snr(snr_db)

    def transmit_cde(self, orig_path, snr_db, sigma_b):
        _, t_path = tempfile.mkstemp()
        _, y_path = tempfile.mkstemp()

        gcx = get_snr(snr_db)
        sigma = gcx * sigma_b
        print(
            f'{self.directory}/s2t -sfile {orig_path} -k 48 -n 48 -Gfile codes/96.3.963/G -smn 1 -tfile {t_path} ')
        print(f'{self.directory}/t2y -tfile {t_path} -yfile {y_path} -gcx {gcx} -seed {np.random.randint(0, 2**31)} -n 96 -sigma {sigma} ')

        os.system(
            f'./s2t -sfile {orig_path} -k 48 -n 48 -Gfile codes/96.3.963/G -smn 1 -tfile {t_path} ')

        os.system(
            f'{self.directory}/t2y -tfile {t_path} -yfile {y_path} -gcx {gcx} -seed {np.random.randint(0, 2**31)} -n 96 -sigma {sigma} ')
        t = np.loadtxt(t_path)
        y = np.loadtxt(y_path)
        os.remove(t_path)
        os.remove(y_path)
        return t, y


if __name__ == '__main__':
    L = LDPCGenerator()
    t, x, y = L(0, 0)

    print(t)
    print(x)
    print(y)
