import numpy as np
from scipy.fftpack import dct
import ipdb
import json
import codecs
from collections import namedtuple


test_data_path = "/home/arnold/matlab/l1c/test/test_data"


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
    print(E_mat.shape)

    Amat = E_mat.dot(dct_mat.T)
    Atmat = Amat.T

    return Amat, Atmat, pix_idx


def f_fun(xu, b, Amat, epsilon, tau):
    n = int(len(xu)/2)

    x = xu[0:n]
    u = xu[n:]

    r = Amat.dot(x) - b
    fu1 = x-u
    fu2 = -x - u
    fe = (1.0/2.0) * (r.T.dot(r) - epsilon**2)
    # ipdb.set_trace()
    f = np.sum(u) - (1.0/tau)*(np.sum(np.log(-fu1)) +
                               np.sum(np.log(-fu2)) + np.log(-fe))

    return r, f, fe, fu1, fu2


def f_fun_fact(b, Amat, epsilon, tau):
    def fun(xu):
        _, fx, _, _, _, = f_fun(xu, b, Amat, epsilon, tau)
        return fx

    return fun


def gradf(Amat, r, fu1, fu2, fe, tau):
    atr = (Amat.T).dot(r)

    ntau_gradx = 1/fu1 - 1/fu2 + (1/fe) * atr
    ntau_gradu = -tau - 1/fu1 - 1/fu2
    gradf = -(1/tau)*np.vstack((ntau_gradx, ntau_gradu))

    return gradf, ntau_gradu, ntau_gradx


def Hess_gradf(fe, fu1, fu2, r, A):

    # gf, ntgu, ntgz = gradf(A, r, fu1, fu2, fe, tau)
    atr = (A.T).dot(r)

    sig11 = (1.0/fu1**2) + (1.0/fu2**2)
    sig12 = -(1.0/fu1**2) + (1.0/fu2**2)

    H11 = np.diag(sig11[:, 0]) - (1.0/fe) * A.T.dot(A) + (1.0/fe**2) * atr.dot(atr.T)

    H1 = np.hstack((H11, np.diag(sig12[:, 0])))
    H2 = np.hstack((np.diag(sig12[:, 0]), np.diag(sig11[:, 0])))
    H = np.vstack((H1, H2))

    return H, H1, H2


def grad_est(x0, h, fun):
    """
    Compute an estimate of the gradient of the functional fun via finite
    difference.
    """
    N = len(x0)

    dx = np.zeros((N, 1))
    ek = np.zeros((N, 1))
    fx0 = fun(x0)

    for k in range(0, N):
        ek = ek*0
        ek[k] = 1.0
        dx[k] = (fun(x0 + ek*h) - fun(x0 - ek*h)) / (2.0*h)

    return dx


def Hess_est(x0, h, fun):
    """
    Compute an estimate of the Hessian of the functional fun via finite
    difference.
    For a function of N variables, we have

    d^2f       f(xj+h, xk+h) - f(xj, xk+h) - f(xj+h, xk) + f(xj, xk)
    -------- =------------------------------------------------------ = H[j,k]
    dxj dxk                   h^2

    where for brevity, the variables with indeces not equal j,k are not shown
    and stay constant.
    """

    N = len(x0)

    H = np.zeros((N, N))
    ek_row = np.zeros((N, 1))
    ek_col = np.zeros((N, 1))
    # fx0 = fun(x0)

    for jcol in range(0, N):
        for krow in range(0, N):
            ek_row = ek_row*0
            ek_col = ek_col*0
            ek_row[krow] = 1.0
            ek_col[jcol] = 1.0
            f12p = fun(x0 + ek_col*h + ek_row*h)
            f12m = fun(x0 - ek_col*h - ek_row*h)
            f1 = fun(x0 + ek_col*h - ek_row*h)
            f2 = fun(x0 - ek_col*h + ek_row*h)

            H[jcol, krow] = ( (f12p+f12m) - (f1 + f2)) / (4.0*h*h)

    return H


np.random.seed(0)
N = 16
tau = 10
mu = 10
lbtol = 1e-3
epsilon = .1
h = 5e-6
A, At, pix_idx = build_Amats(N)
M = A.shape[0]
xi = np.random.randn(N, 1)

b = A.dot(xi) + np.random.randn(M, 1) * 0.015
x = At.dot(b) + np.random.randn(N, 1) * 0.0001


u = (0.95)*np.abs(x) + (0.10)*np.max(np.abs(x))

xu = np.vstack((x, u))
r, f, fe, fu1, fu2 = f_fun(xu, b, A, epsilon, tau)

print(f)
print(fe)


f_fun_handle = f_fun_fact(b, A, epsilon, tau)

gf_true = gradf(A, r, fu1, fu2, fe, tau)[0]

gf_est = grad_est(xu, h, f_fun_handle)


print("gf_true    gf_est   error")
print(np.hstack((gf_true,  gf_est, gf_true-gf_est)))

H_est = Hess_est(xu, h, f_fun_handle)*tau
H_true, H1, H2 = Hess_gradf(fe, fu1, fu2, r, A)

np.set_printoptions(threshold=10000, linewidth=200, precision=4)
print(H_est - H_true)

Herr = H_est - H_true

