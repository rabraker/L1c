#!/usr/bin/env python3
import numpy as np
import json
import codecs
import os
from numpy.random import seed, rand
from scipy.fftpack import dct
import L1cTestDataUtils as TDU


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


def build_dct_rand_test_data(fname, pix_idx, Nx):

    seed(0)

    eta_vec = rand(Nx)

    Adct = Adct_factory(pix_idx)
    Atdct = Atdct_factory(pix_idx, Nx)

    # Another small example with randomly generated data.

    # y = E*M*x
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    # --------------------- M^TE^T EM ----------------------
    # data = struct('x0', x(:)', 'x1', x_AtA(:)', 'pix_idx', pix_idx(:)'-1)
    data = {'x_in': eta_vec.tolist(),
            'y_in': y_vec.tolist(),
            'EMx': EMx.tolist(),
            'pix_idx': pix_idx.tolist(),
            'MtEty': MtEty.tolist(),
            'MtEt_EMx': MtEt_EMx.tolist()}

    TDU.save_json(data, fname)


srcdir = os.getenv("srcdir")
if srcdir is None:
    srcdir = "."

data_dir = srcdir+"/test_data"

fname = data_dir+"/dct_small.json"


Nx = 50
pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])


# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.
build_dct_rand_test_data(fname, pix_idx, Nx)
pix_idx = np.arange(50)
fname = data_dir+"/dct_small_pure_dct.json"
build_dct_rand_test_data(fname, pix_idx, Nx)
