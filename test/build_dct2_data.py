#!/usr/bin/env python3
import codecs
import json
import os

import L1cTestDataUtils as TDU
import numpy as np
from numpy.random import rand, seed
from scipy.fftpack import dct


def Mx_fun(x, mrow, mcol):
    x = x.reshape((mrow, mcol))
    y = dct(x, axis=0, type=3, norm="ortho")
    y = dct(y, axis=1, type=3, norm="ortho")
    y = y.flatten()
    return y


def Mty_fun(z, mrow, mcol):
    z = z.reshape((mrow, mcol))
    x = dct(z, axis=0, type=2, norm="ortho")
    x = dct(x, axis=1, type=2, norm="ortho")
    return x


def Adct_factory(pix_idx, mrow, mcol):
    def Adct(x):
        x = x.reshape((mrow, mcol))
        y = dct(x, axis=0, type=3, norm="ortho")
        y = dct(y, axis=1, type=3, norm="ortho")
        y = y.flatten()[pix_idx]
        return y

    return Adct


def Atdct_factory(pix_idx, mrow, mcol):
    def Atdct(y):
        x = np.zeros(mrow * mcol)
        x[pix_idx] = y
        x = x.reshape((mrow, mcol))
        x = dct(x, axis=0, type=2, norm="ortho")
        x = dct(x, axis=1, type=2, norm="ortho")
        return x

    return Atdct


def build_dct_rand_test_data(fname, pix_idx, mrow, mcol):

    seed(0)

    eta_vec = rand(mrow * mcol)
    z_vec = rand(mrow * mcol)
    Adct = Adct_factory(pix_idx, mrow, mcol)
    Atdct = Atdct_factory(pix_idx, mrow, mcol)

    # Another small example with randomly generated data.

    # y = E*M*x
    y_vec = dct(eta_vec, norm="ortho", type=3)[pix_idx]
    Mx = Mx_fun(eta_vec, mrow, mcol)
    Mty = Mty_fun(z_vec, mrow, mcol)
    Ex = eta_vec[pix_idx]
    Ety = np.zeros_like(eta_vec)
    Ety[pix_idx] = y_vec

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    data = {
        "x_in": eta_vec,
        "y_in": y_vec,
        "z_in": z_vec,
        "EMx": EMx,
        "Mx": Mx,
        "Mty": Mty,
        "Ex": Ex,
        "Ety": Ety,
        "pix_idx": pix_idx,
        "MtEty": MtEty,
        "MtEt_EMx": MtEt_EMx,
        "mrow": mrow,
        "mcol": mcol,
    }

    data = TDU.jsonify(data)
    TDU.save_json(data, fname)


data_dir = TDU.data_dir()

# -------------------------------------------------------- #
fname = data_dir + "/dct2_small.json"
mrow = 16
mcol = 16

pix_idx = np.array([0, 2, 10, 15, 20, 25, 30, 35, 40, 45, 49])
build_dct_rand_test_data(fname, pix_idx, mrow, mcol)


# -------------------------------------------------------- #
# Check we still do this right for a tall, skinny matrix
fname = data_dir + "/dct2_small_tall.json"
mrow = 18
mcol = 16

# pix_idx = np.array([0, 2, 10, 15, 20, 25, 30, 35, 40, 45, (mrow-5)*mcol])
pix_idx = np.arange(mrow * mcol)
build_dct_rand_test_data(fname, pix_idx, mrow, mcol)

# -------------------------------------------------------- #
# Check we still do this right for a wide, short matrix
fname = data_dir + "/dct2_small_wide.json"
mrow = 16
mcol = 18

# pix_idx = np.array([0, 2, 10, 15, 20, 25, 30, 35, 40, 45, (mrow-5)*mcol])
pix_idx = np.arange(mrow * mcol)
build_dct_rand_test_data(fname, pix_idx, mrow, mcol)


# -------------------------------------------------------- #
# This checks that we get the right result for a pure dct, and in
# particular, that we got the scaling at the first element correct.

mrow = 16
mcol = 16

pix_idx = np.arange(mrow * mcol)
fname = data_dir + "/dct2_small_pure_dct.json"
build_dct_rand_test_data(fname, pix_idx, mrow, mcol)
