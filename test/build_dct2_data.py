#!/usr/bin/env python3
import numpy as np
import json
import codecs
import os
from numpy.random import seed, rand
from scipy.fftpack import dct


def save_json(data, fname):
    json.dump(data, codecs.open(fname, 'w', encoding='utf-8'),
              separators=(',', ':'), sort_keys=True, indent=4)


def Adct_factory(pix_idx, N, M):
    def Adct(x):
        x = x.reshape((N, M))
        y = dct(x, axis=0, type=3, norm='ortho')
        y = dct(y, axis=1, type=3, norm='ortho')
        y = y.flatten()[pix_idx]
        return y
    return Adct


def Atdct_factory(pix_idx, N, M):
    def Atdct(y):
        x = np.zeros(N*M)
        x[pix_idx] = y
        x = x.reshape((N, M))
        x = dct(x, axis=0, type=2, norm='ortho')
        x = dct(x, axis=1, type=2, norm='ortho')
        return x

    return Atdct


def build_dct_rand_test_data(fname, pix_idx, N, M):

    seed(0)

    eta_vec = rand(N*M)

    Adct = Adct_factory(pix_idx, N, M)
    Atdct = Atdct_factory(pix_idx, N, M)

    # Another small example with randomly generated data.

    # y = E*M*x
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    data = {'x_in': eta_vec.flatten().tolist(),
            'y_in': y_vec.flatten().tolist(),
            'EMx': EMx.flatten().tolist(),
            'pix_idx': pix_idx.flatten().tolist(),
            'MtEty': MtEty.flatten().tolist(),
            'MtEt_EMx': MtEt_EMx.flatten().tolist(),
            'N': N,
            'M': M}

    save_json(data, fname)


srcdir = os.getenv("srcdir")
if srcdir is None:
    srcdir = "."

data_dir = srcdir+"/test_data"

fname = data_dir+"/dct2_small.json"


N = 16
M = 16

pix_idx = np.array([0, 2, 10, 15, 20, 25, 30, 35, 40, 45, 49])
build_dct_rand_test_data(fname, pix_idx, N, M)


# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.


pix_idx = np.arange(N*M)
fname = data_dir+"/dct2_small_pure_dct.json"
build_dct_rand_test_data(fname, pix_idx, N, M)
