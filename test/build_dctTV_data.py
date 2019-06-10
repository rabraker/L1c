#!/usr/bin/env python3
import numpy as np
import json
import codecs
import os
from numpy.random import seed, rand
from scipy.fftpack import dct

import L1cTestDataUtils as TDU
import build_dct2_data as dct2
# import build_dct_data as dct1
import build_TV_data as TV

import ipdb


def Wz_factory(mrow, mcol, alp_v, alp_h):
    Dy_mat = TV.DyMatRep(mrow, mcol)
    Dx_mat = TV.DxMatRep(mrow, mcol)

    def Wx_fun(z):
        """
        Wx operator for both Dx and Dv and dct2
        """
        m = mrow * mcol
        u0 = z[0:m]

        x = dct2.Mx_fun(u0, mrow, mcol).flatten()
        if alp_v > 0:
            u1 = z[m:2*m]
            x = x + alp_v * Dy_mat.dot(u1).flatten()
        if alp_h > 0:
            if alp_v > 0:
                u2 = z[2*m:3*m]
            else:
                u2 = z[m:2*m]
            x = x + alp_h * Dx_mat.dot(u2).flatten()

        return x

    return Wx_fun


def Wtx_factory(mrow, mcol, alp_v, alp_h):
    Dy_mat = TV.DyMatRep(mrow, mcol)
    Dx_mat = TV.DxMatRep(mrow, mcol)

    def Wtx_fun(x):
        """
        Wx operator for both Dx and Dv and dct2
        """
        y1 = dct2.Mty_fun(x, mrow, mcol).flatten()
        y2 = alp_v * Dy_mat.T.dot(x).flatten()
        y3 = alp_h * Dx_mat.T.dot(x).flatten()

        if alp_v > 0 and alp_h > 0:
            return np.hstack((y1, y2, y3))

        if alp_v <= 0 and alp_h > 0:
            return np.hstack((y1, y3))

        if alp_v > 0 and alp_h <= 0:
            return np.hstack((y1, y2))

        return y1

    return Wtx_fun


def Az_factory(pix_idx, Wz_fun):
    def Az_fun(z):
        x = Wz_fun(z)
        y = x[pix_idx]
        return y

    return Az_fun


def Aty_factory(pix_idx, Wtx_fun, m):
    def Aty_fun(y):
        x = np.zeros(m)
        x[pix_idx] = y
        z = Wtx_fun(x)
        return z

    return Aty_fun


def build_dct2_TV_vh(fname, pix_idx, mrow, mcol, alp_v, alp_h):
    """
    For simplicity, take
    W = [M, Dx], and W^T = [M^T;
                            Dx^T];
    And E is n by m,
        W is m by p, p=2*m

    Then,
    x = W * z
    y = E * W * z = E * x

    so z is in R^p, x in R^m
    """

    seed(0)

    n = len(pix_idx)

    m = mrow*mcol
    p = m
    if alp_v > 0:
        p = p + m
    if alp_h > 0:
        p = p + m

    x_vec = rand(mrow*mcol)
    z_vec = rand(p)

    Wz_fun = Wz_factory(mrow, mcol, alp_v, alp_h)
    Wtx_fun = Wtx_factory(mrow, mcol, alp_v, alp_h)

    A = Az_factory(pix_idx, Wz_fun)
    At = Aty_factory(pix_idx, Wtx_fun, mrow*mcol)

    # Another small example with randomly generated data.

    # y = E*M*x
    y_vec = rand(n)
    Mx = Wz_fun(z_vec)
    Mty = Wtx_fun(x_vec)
    Ex = x_vec[pix_idx]
    Ety = np.zeros_like(x_vec)
    Ety[pix_idx] = y_vec

    EMx = A(z_vec)
    MtEty = At(y_vec)
    MtEt_EMx = At(A(z_vec))

    data = {'x_in': x_vec,
            'y_in': y_vec,
            'z_in': z_vec,
            'EMx': EMx,
            'Mx': Mx,
            'Mty': Mty,
            'Ex': Ex,
            'Ety': Ety,
            'pix_idx': pix_idx,
            'MtEty': MtEty,
            'MtEt_EMx': MtEt_EMx,
            'alp_v': alp_v,
            'alp_h': alp_h,
            'mrow': mrow,
            'mcol': mcol,
            'p': p}

    data = TDU.jsonify(data)
    TDU.save_json(data, fname)


if __name__ == "__main__":
    srcdir = os.getenv("srcdir")
    if srcdir is None:
        srcdir = "."

    data_dir = srcdir+"/test_data"

    # -------------------------------------------------------- #
    fname = data_dir+"/dct2_tv_square.json"
    mrow = 16
    mcol = 16

    pix_idx = np.array([0, 2, 10, 15, 20, 25, 30, 35, 40, 45, 49])
    alp_v = 0
    alp_h = 0
    build_dct2_TV_vh(fname, pix_idx, mrow, mcol, alp_v, alp_h)

    # -------------------------------------------------------- #
    fname = data_dir+"/dct2_tv_vh_square.json"
    alp_v = 1.5
    alp_h = 2.0
    build_dct2_TV_vh(fname, pix_idx, mrow, mcol, alp_v, alp_h)

    # -------------------------------------------------------- #
    fname = data_dir+"/dct2_tv_v_square.json"
    alp_v = 1.5
    alp_h = 0.0
    build_dct2_TV_vh(fname, pix_idx, mrow, mcol, alp_v, alp_h)

    # -------------------------------------------------------- #
    fname = data_dir+"/dct2_tv_h_square.json"
    alp_v = 0.0
    alp_h = 2.0
    build_dct2_TV_vh(fname, pix_idx, mrow, mcol, alp_v, alp_h)
