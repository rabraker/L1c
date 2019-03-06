#!/usr/bin/env python3
import numpy as np
import numpy.testing as npt
import TestDataUtils as TDU
import build_TV_data as TV
import ipdb


def shrink1(x, gamma):

    zr = np.zeros_like(x)
    y = np.sign(x) * np.max((np.abs(x) - gamma, zr), axis=0)
    return y


def norm2_err(x, y):

    err = x - y

    return err.T.dot(err)


def breg_rhs(n, m, lam, mu, f, dx, bx, dy, by):
    DxMat = TV.DxMatRep(n, m)
    DyMat = TV.DyMatRep(n, m)

    dxx = lam * DxMat.T.dot(dx-bx)
    dyy = lam * DyMat.T.dot(dy-by)

    return mu*f + dxx + dyy

def breg_hess_eval(n, m, mu, lam, x):
    DxMat = TV.DxMatRep(n, m)
    DyMat = TV.DyMatRep(n, m)

    dxx = DxMat.T.dot(DxMat).dot(x)
    dyy = DyMat.T.dot(DyMat).dot(x)

    # ipdb.set_trace()
    return mu - lam*(dxx + dyy)


def build_data(fname):
    np.random.seed(0)
    n = 4
    m = 4
    N = n*m

    gamma = .53
    lam = 1
    mu = 1

    x = np.random.rand(N, 1)
    y = np.random.rand(N, 1)

    dx = np.random.rand(N, 1)
    bx = np.random.rand(N, 1)
    dy = np.random.rand(N, 1)
    by = np.random.rand(N, 1)
    f = np.random.rand(N, 1)

    rhs = breg_rhs(n, m, lam, mu, f, dx, bx, dy, by)

    Hess_x = breg_hess_eval(n, m, mu, lam, x)

    x_shrunk = shrink1(x, gamma)
    z = -x + y

    nrm2 = norm2_err(x, y)

    data = {'N': N,
            'n': n,
            'm': m,
            'gamma': gamma,
            'x': x,
            'y': y,
            'z': z,
            'dx': dx,
            'bx': bx,
            'dy': dy,
            'by': by,
            'f': f,
            'lam': lam,
            'mu': mu,
            'rhs': rhs,
            'Hessx': Hess_x,
            'x_shrunk': x_shrunk,
            'nrm2': nrm2}

    data = TDU.jsonify(data)

    TDU.save_json(data, fname)


if __name__ == '__main__':
    import os

    srcdir = os.getenv("srcdir")
    if srcdir is None:
        srcdir = "."

    data_dir = srcdir+"/test_data"

    fname = data_dir+"/bregman.json"

    build_data(fname)
