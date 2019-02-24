#!/usr/bin/env python3
import os
import sys
import numpy as np
from scipy.fftpack import dct
import ipdb
import json
import codecs


def jsonify(data):

    for key, value in data.items():
        if type(value) is np.ndarray:
            if len(value.flatten()) == 1:
                data[key] = value.flatten()[0]
            else:
                data[key] = value.flatten().tolist()

    return data


def save_json(data, file_path):

    json.dump(data, codecs.open(file_path, 'w'),
              separators=(',', ':'), sort_keys=True, indent=4)


def build_dct_mat(N):
    dct_mat = np.zeros((N, N))
    I = np.eye(N)

    for k in range(0, N):
        dct_mat[:, k] = dct(I[:, k], type=2, norm='ortho')

    return dct_mat


def build_Amats(N):

    dct_mat = build_dct_mat(N)

    E_mat = np.eye(N)
    row_select = np.random.rand(N)

    # The comma is important! Otherwise we get tuple and everything
    # below is fubar.
    pix_idx, = np.where(row_select > 0.7)

    E_mat = E_mat[pix_idx, :]

    Amat = E_mat.dot(dct_mat.T)
    Atmat = Amat.T

    return Amat, Atmat, pix_idx


def find_max_step(dx, du, Adx, fu1, fu2, r, epsilon):
    idx_fu1, = np.where((dx-du)[:, 0] > 0)
    idx_fu2, = np.where((-dx-du)[:, 0] > 0)

    aqe = Adx.T.dot(Adx)[0][0]
    bqe = 2*r.T.dot(Adx)[0][0]
    cqe = (r.T.dot(r) - epsilon**2)[0][0]

    smax_f1 = np.min(-fu1[idx_fu1] / (dx[idx_fu1] - du[idx_fu1]))
    smax_f2 = np.min(-fu2[idx_fu2] / (-dx[idx_fu2] - du[idx_fu2]))
    smax_quad = (-bqe + np.sqrt(bqe**2 - 4*aqe*cqe)) / (2*aqe)

    smax = np.min([1, smax_f1, smax_f2, smax_quad])

    return smax * 0.99


def newton_init(x0, lbtol, mu):
    N = x0.shape[0]
    u = (0.95)*np.abs(x0) + (0.10)*np.max(np.abs(x0))
    tmp = (2.0*N+1.0)/np.sum(np.abs(x0))
    tau = np.max((tmp, 1))
    lbiter = np.ceil((np.log(2.0*N+1.0)-np.log(lbtol)-np.log(tau))/np.log(mu))

    return u, tau, lbiter


def f_fun(xu, b, Amat, epsilon, tau):
    n = int(len(xu)/2)

    x = xu[0:n]
    u = xu[n:]

    r = Amat.dot(x) - b
    fu1 = x-u
    fu2 = -x - u
    fe = (1.0/2.0) * (r.T.dot(r) - epsilon**2)
    f = np.sum(u) - (1.0/tau)*(np.sum(np.log(-fu1)) +
                               np.sum(np.log(-fu2)) + np.log(-fe))

    return r, f, fe, fu1, fu2


def f_fun_fact(b, Amat, epsilon, tau):
    def fun(xu):
        _, fx, _, _, _, = f_fun(xu, b, Amat, epsilon, tau)
        return fx

    return fun


def gradf(Amat, atr, fu1, fu2, fe, tau):
    # atr = (Amat.T).dot(r)

    ntau_gradx = 1/fu1 - 1/fu2 + (1/fe) * atr
    ntau_gradu = -tau - 1/fu1 - 1/fu2
    gradf = -(1/tau)*np.vstack((ntau_gradx, ntau_gradu))

    return gradf, ntau_gradu, ntau_gradx


def H11p_fun(sigx, fe, A, atr, z):
    # Need the slash here, python will not continue the line on its own.
    h11p_z = sigx*z + (-1.0/fe) * A.T.dot(A.dot(z)) + \
             (1/fe**2) * (atr.T.dot(z))*atr

    return h11p_z


def Hess_gradf(fe, fu1, fu2, atr, A):

    gf, ntgu, ntgz = gradf(A, atr, fu1, fu2, fe, tau)
    # atr = (A.T).dot(r)

    sig11 = (1.0/fu1**2) + (1.0/fu2**2)
    sig12 = -(1.0/fu1**2) + (1.0/fu2**2)

    H11 = np.diag(sig11[:, 0]) - (1.0/fe) * A.T.dot(A) + (1.0/fe**2) * atr.dot(atr.T)

    H1 = np.hstack((H11, np.diag(sig12[:, 0])))
    H2 = np.hstack((np.diag(sig12[:, 0]), np.diag(sig11[:, 0])))
    H = np.vstack((H1, H2))

    sigx = sig11 - (sig12**2)/sig11
    w1p = ntgz - (sig12/sig11)*ntgu

    H11_prime = np.diag(sigx[:, 0]) - (1.0/fe) * A.T.dot(A) + \
        (1.0/fe**2) * atr.dot(atr.T)

    dxdu = np.linalg.solve(H, np.vstack((ntgz, ntgu)))

    dx = np.linalg.solve(H11_prime, w1p)
    du = (1/sig11)*ntgu - (sig12/sig11)*dx
    np.testing.assert_array_almost_equal(dx, dxdu[0:N])
    np.testing.assert_array_almost_equal(du, dxdu[N:])

    return gf, ntgu, ntgz, atr, sig11, sig12, sigx, w1p, dxdu


srcdir = os.getenv("srcdir")
if srcdir is None:
    srcdir = "."

test_data_path = srcdir+"/test_data"

l1qc_data_path = test_data_path+"/l1qc_data.json"

np.random.seed(0)
N = 16
tau = 10
mu = 10
lbtol = 1e-3
epsilon = .1

A, At, pix_idx = build_Amats(N)
M = A.shape[0]

xx = np.random.randn(N, 1)
b = A.dot(xx) + np.random.randn(M, 1) * 0.015
x = At.dot(b) + np.random.randn(N, 1) * 0.0001

u, tau_exp, lbiter = newton_init(x, lbtol, mu)

xu = np.vstack((x, u))
r, f, fe, fu1, fu2 = f_fun(xu, b, A, epsilon, tau)

if np.isnan(f) or np.isnan(fe):
    print("Randomly generated test data is infeasible. ")
    sys.exit(1)


dx_rand1 = np.random.randn(N, 1)
du_rand1 = np.random.randn(N, 1)
Adx_rand1 = A.dot(dx_rand1)
smax = find_max_step(dx_rand1, du_rand1, Adx_rand1, fu1, fu2, r, epsilon)

atr = A.T.dot(r)
gf, ntau_gradu, ntau_gradx = gradf(A, atr, fu1, fu2, fe, tau)
gf, ntgu, ntgx, atr, sig11, sig12, sigx, w1p, dxdu = Hess_gradf(fe, fu1,
                                                                fu2, atr, A)

sigx_rand = np.random.randn(N, 1)
z_rand = np.random.randn(N, 1)
r_rand = np.random.randn(M, 1)
at_rrand = A.T.dot(r)
fe_rand = -np.abs(np.random.rand(1)[0])
h11p_z = H11p_fun(sigx_rand, fe_rand, A, at_rrand, z_rand)

l1qc_data = {'N': N,
             'M': M,
             'tau0': tau,
             'lbtol': lbtol,
             'mu': mu,
             'epsilon': epsilon,
             'cgtol': 1e-9,
             'cgmaxiter': 500,
             'pix_idx': pix_idx,
             'A': A,
             'b': b,
             'x': x,
             # newton init
             'tau': tau_exp,
             'lbiter': lbiter,
             'u': u,
             'r': r,
             'f': f,
             'fe': fe,
             'fu1': fu1,
             'fu2': fu2,
             # for find_max_step
             'Adx_rand1': Adx_rand1,
             'dx_rand1': dx_rand1,
             'du_rand1': du_rand1,
             'smax': smax,
             # for h11p
             'sigx_rand': sigx_rand,
             'z_rand': z_rand,
             'r_rand': r_rand,
             'at_rrand': at_rrand,
             'fe_rand': fe_rand,
             'h11p_z': h11p_z,
             # for decent data
             'gradf': gf,
             'atr': atr,
             'ntgx': ntgx,
             'ntgu': ntgu,
             'sig11': sig11,
             'sig12': sig12,
             'w1p': w1p,
             'sigx': sigx,
             'dx': dxdu[0:N],
             'du': dxdu[N:]}


l1qc_data_json = jsonify(l1qc_data)

save_json(l1qc_data, l1qc_data_path)
