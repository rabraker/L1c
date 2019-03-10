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

    y_nrm = y.T.dot(y)
    err = x - y
    err_nrm = err.T.dot(err)
    return np.sqrt(err_nrm/y_nrm)


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
    return mu*x + lam*(dxx + dyy)


def breg_hess_solve(n, m, mu, lam, b):
    DxMat = TV.DxMatRep(n, m)
    DyMat = TV.DyMatRep(n, m)

    Dxx = DxMat.T.dot(DxMat)
    Dyy = DyMat.T.dot(DyMat)
    I = np.eye(n*m)

    H = mu*I + lam*(Dxx + Dyy)
    return np.linalg.solve(H, b)


def Hess_diag(N, M, mu, lam):
    import build_TV_data
    DyMat = build_TV_data.DyMatRep(N, M)
    DxMat = build_TV_data.DxMatRep(N, M)

    H = DyMat.T.dot(DyMat) + DxMat.T.dot(DxMat)
    I = np.eye(N*M)
    H = mu*I + lam*H
    H = np.diag(H)
    H = 1/H
    return H


def build_data(fname):
    np.random.seed(0)
    n = 16
    m = 16
    N = n*m


    gamma = .53
    mu = 10
    lam = 2*mu

    hdiag = Hess_diag(n, m, mu, lam)

    x = np.random.rand(N, 1)
    y = np.random.rand(N, 1)

    dx = np.random.rand(N, 1)
    bx = np.random.rand(N, 1)
    dy = np.random.rand(N, 1)
    by = np.random.rand(N, 1)
    f = np.random.rand(N, 1)
    b = np.random.rand(N, 1)
    Hsolve_b = breg_hess_solve(n, m, mu, lam, b)

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
            'b': b,
            'Hsolveb': Hsolve_b,
            'Hessx': Hess_x,
            'H_diag': hdiag,
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
