#!/usr/bin/env python3
import numpy as np
import os
from numpy.random import seed, rand
from scipy.fftpack import dct
import L1cTestDataUtils as TDU
import build_CS20NG_example_data as CS20NG


def Adct_factory(pix_idx):
    def Adct(x):
        y = dct(x, type=3, norm='ortho')
        y = y[pix_idx]
        return y
    return Adct


def Atdct_factory(pix_idx, m):
    def Atdct(y):
        x = np.zeros(m)
        x[pix_idx] = y
        x = dct(x, type=2, norm='ortho')
        return x

    return Atdct


def build_dct_test_data(fname, pix_idx, m, eta_vec=None):
    """
    Build a small example with random test data for the eta vector.
    y = EM * eta
    x = M*eta
    """

    seed(0)
    if eta_vec is None:
        eta_vec = rand(m)

    z_vec = rand(m)

    Adct = Adct_factory(pix_idx)
    Atdct = Atdct_factory(pix_idx, m)

    # Another small example with randomly generated data.

    # y = E*M*eta
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]
    Mx = dct(eta_vec, type=3, norm='ortho')
    Mty = dct(z_vec, type=2, norm='ortho')
    Ex = eta_vec[pix_idx]
    Ety = np.zeros_like(eta_vec)
    Ety[pix_idx] = y_vec

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    data = {'x_in': eta_vec,
            'y_in': y_vec,
            'z_in': z_vec,
            'EMx': EMx,
            'Mx': Mx,
            'Mty': Mty,
            'Ex': Ex,
            'Ety': Ety,
            'pix_idx': pix_idx,
            'MtEty': MtEty,
            'MtEt_EMx': MtEt_EMx}

    data = TDU.jsonify(data)
    TDU.save_json(data, fname)


def cs20ng_example(npix):
    """
    Build a small example derived from the simulated CS20NG sample grating,
    with a mu-path mask applied.
    """
    m = npix ** 2
    perc = 0.15
    mu_len = 25
    img = CS20NG.make_CS20NG(npix)
    pix_idx, _ = CS20NG.mu_path_mask(mu_len, npix, perc)

    x = img.flatten()
    eta = dct(x, norm='ortho', type=2)

    return eta,  pix_idx, m


data_dir = TDU.data_dir()

# ------------- A small example, with sub-sampling ------------
fname = data_dir+"/dct_small.json"
m = 50
pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])
build_dct_test_data(fname, pix_idx, m)

# ------------- A small example, without sub-sampling ------------
# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.
pix_idx = np.arange(50)
fname = data_dir+"/dct_small_pure_dct.json"
build_dct_test_data(fname, pix_idx, m)


# -------- A large set of test data -------------
fname = data_dir+"/dct_large.json"
npix = 256
eta, pix_idx, m = cs20ng_example(npix)
build_dct_test_data(fname, pix_idx, m, eta_vec=eta)
