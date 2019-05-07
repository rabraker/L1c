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


def dct_test_data(m, eta_vec, z_vec, pix_idx):
    Adct = Adct_factory(pix_idx)
    Atdct = Atdct_factory(pix_idx, m)

    # Another small example with randomly generated data.

    # y = E*M*eta
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]

    Mx = dct(eta_vec, type=3, norm='ortho')
    Mty = dct(z_vec, type=2, norm='ortho')
    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    return y_vec, Mx, Mty, EMx, MtEty, MtEt_EMx


def save_dct_test_data(fname, eta, y, z, Mx, Mty, EMx, MtEty, MtEt_EMx, pix_idx):
    data = {'x_in': eta,
            'y_in': y,
            'z_in': z,
            'EMx': EMx,
            'Mx': Mx,
            'Mty': Mty,
            'pix_idx': pix_idx,
            'MtEty': MtEty,
            'MtEt_EMx': MtEt_EMx}
    data = TDU.jsonify(data)
    TDU.save_json(data, fname)


def build_dct_rand_test_data(fname, pix_idx, m):
    """
    Build a small example with random test data for the eta vector.
    y = EM * eta
    x = M*eta
    """

    seed(0)
    eta_vec = rand(m)
    z = rand(m)
    y_vec, Mx, Mty, EMx, MtEty, MtEt_EMx = dct_test_data(m, eta_vec, z, pix_idx)

    save_dct_test_data(fname, eta_vec, y_vec, z, Mx, Mty, EMx, MtEty, MtEt_EMx, pix_idx)


def build_dct_large(fname, npix):
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
    z = rand(m);

    y_vec, Mx, Mty, EMx, MtEty, MtEt_EMx = dct_test_data(m, eta, z, pix_idx)

    save_dct_test_data(fname, eta, y_vec, z, Mx, Mty, EMx, MtEty, MtEt_EMx, pix_idx)


srcdir = os.getenv("srcdir")
if srcdir is None:
    srcdir = "."

data_dir = srcdir+"/test_data"

# ------------- A small example, with sub-sampling ------------
fname = data_dir+"/dct_small.json"
m = 50
pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])

# ------------- A small example, without sub-sampling ------------
# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.
build_dct_rand_test_data(fname, pix_idx, m)
pix_idx = np.arange(50)
fname = data_dir+"/dct_small_pure_dct.json"
build_dct_rand_test_data(fname, pix_idx, m)


# -------- A large set of test data -------------
fname = data_dir+"/dct_large.json"
npix = 256
build_dct_large(fname, npix)
