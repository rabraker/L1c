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


def Atdct_factory(pix_idx, N):
    def Atdct(y):
        x = np.zeros(N)
        x[pix_idx] = y
        x = dct(x, type=2, norm='ortho')
        return x

    return Atdct


def dct_test_data(Nx, eta_vec, pix_idx):
    Adct = Adct_factory(pix_idx)
    Atdct = Atdct_factory(pix_idx, Nx)

    # Another small example with randomly generated data.

    # y = E*M*eta
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    return y_vec, EMx, MtEty, MtEt_EMx


def save_dct_test_data(fname, eta, y, EMx, MtEty, MtEt_EMx, pix_idx):
    data = {'x_in': eta.flatten().tolist(),
            'y_in': y.flatten().tolist(),
            'EMx': EMx.flatten().tolist(),
            'pix_idx': pix_idx.flatten().tolist(),
            'MtEty': MtEty.flatten().tolist(),
            'MtEt_EMx': MtEt_EMx.flatten().tolist()}

    TDU.save_json(data, fname)


def build_dct_rand_test_data(fname, pix_idx, Nx):
    """
    Build a small example with random test data for the eta vector.
    y = EM * eta
    x = M*eta
    """

    seed(0)
    eta_vec = rand(Nx)
    y_vec, EMx, MtEty, MtEt_EMx = dct_test_data(Nx, eta_vec, pix_idx)

    save_dct_test_data(fname, eta_vec, y_vec, EMx, MtEty, MtEt_EMx, pix_idx)


def build_dct_large(fname, npix):
    """
    Build a small example derived from the simulated CS20NG sample grating,
    with a mu-path mask applied.
    """
    Nx = npix ** 2
    perc = 0.15
    mu_len = 25
    img = CS20NG.make_CS20NG(npix)
    pix_idx, _ = CS20NG.mu_path_mask(mu_len, npix, perc)

    x = img.flatten()
    eta = dct(x, norm='ortho', type=2)

    y_vec, EMx, MtEty, MtEt_EMx = dct_test_data(Nx, eta, pix_idx)

    save_dct_test_data(fname, eta, y_vec, EMx, MtEty, MtEt_EMx, pix_idx)


srcdir = os.getenv("srcdir")
if srcdir is None:
    srcdir = "."

data_dir = srcdir+"/test_data"

# ------------- A small example, with sub-sampling ------------
fname = data_dir+"/dct_small.json"
Nx = 50
pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])

# ------------- A small example, without sub-sampling ------------
# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.
build_dct_rand_test_data(fname, pix_idx, Nx)
pix_idx = np.arange(50)
fname = data_dir+"/dct_small_pure_dct.json"
build_dct_rand_test_data(fname, pix_idx, Nx)


# -------- A large set of test data -------------
fname = data_dir+"/dct_large.json"
npix = 256
build_dct_large(fname, npix)
