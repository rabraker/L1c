#!/usr/bin/env python3

import os
import numpy as np
import numpy.testing as npt
import L1cTestDataUtils as TDU
from scipy.fftpack import dct


def Adct_factory(pix_idx):
    def Adct(x):
        # y = dct(x, type=3, norm='ortho')
        y = x[pix_idx]
        return y
    return Adct


def Atdct_factory(pix_idx, m_):
    def Atdct(y):
        # x = dct(y, type=2, norm='ortho')
        x = np.zeros(m_)
        x[pix_idx] = y
        return x

    return Atdct


def Udct_factory():
    def Adct(x):
        y = dct(x, type=3, norm='ortho')
        return y
    return Adct


def Utdct_factory():
    def Atdct(y):
        x = dct(y, type=2, norm='ortho')
        return x

    return Atdct


def nesta_project(x, g, b, L_mu, A, At, epsilon):
    """
    Implements the closed form solution to

    yk = argmin_{x\\in Q_p} (Lmu/2) || xk - x||_2^2 + <\\Nabla f_mu(xk), x-xk> (3.3)
    where Q_p = ||b - Ayk||_2^2 <= \\epsilon

    i.e., eq (3.5) and (3.6) and (3.7)

    Note that because A is an othogonal projector, then AA'=I (but not A'A).
    Thus, in (3.6), A'AA'b = A'b.
    """

    q = x - (1.0 / L_mu) * g
    Aq = A(q)
    AtAq = At(Aq)
    Atb = At(b)

    nrm_err = np.linalg.norm(b - Aq)
    a0 = L_mu*(nrm_err/epsilon - 1.0)
    lam = np.max((0.0, a0))
    a1 = lam/(L_mu + lam)

    vk = (lam / L_mu) * (1.0 - a1) * Atb + q - a1*AtAq

    return vk, lam


def nesta_f_eval(W, Wt, xk, mu):
    Wtxk = Ut(xk)

    u = Wtxk / np.fmax(mu, np.abs(Wtxk))

    gradf = W(u)
    fx = u.dot(Wtxk) - 0.5 * mu * np.linalg.norm(u)**2

    return fx, gradf

if __name__ == "__main__":

    data_dir = TDU.data_dir()

    pix_idx = np.array([0, 3, 4, 7, 8])
    M = 10
    N = len(pix_idx)

    L = 1.0
    sigma = 0.1
    mu = 0.1
    Lmu = L / mu
    np.random.seed(0)
    xk = np.random.randn(M)
    b = np.random.randn(N)
    g = np.random.randn(M)

    q = xk - (1.0/Lmu)*g

    A = Adct_factory(pix_idx)
    At = Atdct_factory(pix_idx, M)
    U = Udct_factory()
    Ut = Utdct_factory()

    Atb = At(b)


    yk, lam = nesta_project(xk, g, b, Lmu, A, At,  sigma)

    fx, gradf = nesta_f_eval(U, Ut, xk, mu)

    lhs = yk + (lam/Lmu)*At(A(yk))
    rhs = (lam/Lmu)*Atb + xk - (1.0/Lmu)*g

    # print(np.linalg.norm(b - A(yk)) - sigma)
    npt.assert_almost_equal(np.linalg.norm(b - A(yk)) - sigma, 0)
    npt.assert_array_almost_equal(lhs, rhs)
    # print(lhs-rhs)

    q = At(A(q))
    dat = {'N': N,
           'xk': xk,
           'b': b,
           'g': g,
           'q': q,
           'gradf_exp': gradf,
           'fx_exp': fx,
           'yk_exp': yk,
           'pix_idx': pix_idx,
           'Lmu': Lmu,
           'L': L,
           'mu': mu,
           'sigma': sigma,
           }

    dat = TDU.jsonify(dat)

    TDU.save_json(dat, data_dir+"/nesta_data.json")
