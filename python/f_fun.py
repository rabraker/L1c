import numpy as np
from scipy.fftpack import dct
import ipdb

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
    E_mat = E_mat[row_select > 0.7, :]
    print(E_mat.shape)

    Amat = E_mat.dot(dct_mat.T)
    Atmat = Amat.T

    return Amat, Atmat


def f_fun(xu, b, Amat, epsilon, tau):
    n = int(len(xu)/2)

    x = xu[0:n]
    u = xu[n:]

    r = Amat.dot(x) - b
    fu1 = x-u
    fu2 = -x - u
    fe = (1.0/2.0) * (r.dot(r) - epsilon**2)
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

    ntgz = 1/fu1 - 1/fu2 + (1/fe) * atr
    ntgu = -tau - 1/fu1 - 1/fu2
    gradf = -(1/tau)*np.hstack((ntgz, ntgu))

    return gradf




np.random.seed(0)
N = 15
tau = 10
epsilon = .1
Amat, Atmat = build_Amats(N)

# Amat = Amat*0
# Atmat = Atmat*0
M = Amat.shape[0]
h = 1e-7


xx = np.random.randn(N)
b = Amat.dot(xx) + np.random.randn(M) * 0.015
x = Atmat.dot(b)
u = (0.99)*np.abs(x) + (0.10)*np.max(np.abs(x))


xu = np.hstack((x, u))
r, f, fe, fu1, fu2 = f_fun(xu, b, Amat, epsilon, tau)

f_fun_handle = f_fun_fact(b, Amat, epsilon, tau)



gf_true = gradf(Amat, r, fu1, fu2, fe, tau)




# def grad_est(x0, h, fun):
#     """
#     Compute an estimate of the gradient of the functional fun via finite
#     difference.
#     """
#     N = len(x0)

#     dx = np.zeros(N)
#     ek = np.zeros(N)
#     fx0 = fun(x0)

#     for k in range(0, N):
#         ek = ek*0
#         ek[k] = 1.0

#         dx[k] = (fun(x0 + ek*h) - fun(x0 - ek*h)) / (2.0*h)

#     return dx


# def Hess_est(x0, h, fun):
#     """
#     Compute an estimate of the Hessian of the functional fun via finite
#     difference.
#     For a function of N variables, we have

#     d^2f       f(xj+h, xk+h) - f(xj, xk+h) - f(xj+h, xk) + f(xj, xk)
#     -------- =------------------------------------------------------ = H[j,k]
#     dxj dxk                   h^2

#     where for brevity, the variables with indeces not equal j,k are not shown
#     and stay constant.
#     """

#     N = len(x0)

#     H = np.zeros((N, N))
#     ek_row = np.zeros(N)
#     ek_col = np.zeros(N)
#     # fx0 = fun(x0)

#     for jcol in range(0, N):
#         for krow in range(0, N):
#             ek_row = ek_row*0
#             ek_col = ek_col*0
#             ek_row[krow] = 1.0
#             ek_col[jcol] = 1.0
#             f12p = fun(x0 + ek_col*h + ek_row*h)
#             f12m = fun(x0 - ek_col*h - ek_row*h)
#             f1 = fun(x0 + ek_col*h - ek_row*h)
#             f2 = fun(x0 - ek_col*h + ek_row*h)

#             H[jcol, krow] = ( (f12p+f12m) - (f1 + f2)) / (4.0*h*h)

#     return H


# def Hess(fe, fu1, fu2, r, A):
#     atr = np.array([((A.T).dot(r)).T]).T

#     sig11 = np.diag((1.0/fu1**2) + (1.0/fu2**2))
#     sig12 = np.diag(-(1.0/fu1**2) + (1.0/fu2**2))

#     H11 = sig11 - (1.0/fe) * (A.T).dot(A)
#     + (1.0/fe**2) * atr.dot(atr.T)

#     H1 = np.hstack((H11, sig12))
#     H2 = np.hstack((sig12, sig11))

#     H = np.vstack((H1, H2))

#     return H, H11




# # def f2(x, H):
# #     # H = np.eye(len(x))
# #     return x.dot(H).dot(x)


# def f2_fact(H):
#     def f2_fun(x):
#         return 0.5*x.dot(H).dot(x)

#     return f2_fun
#

# gf_est = grad_est(xu, h, f_fun_handle)

# np.set_printoptions(linewidth=300, precision=3)
# print(gf_est-gf_true)

# H, H11 = Hess(fe, fu1, fu2, r, Amat)
# Hest = tau*Hess_est(xu, h, f_fun_handle)


# print(H11)
# print("----------------------\n")
# print(Hest[0:N, 0:N])

# # print(H[0:N, 0:N] -Hest[0:N, 0:N])

# print(np.max(np.abs(H-Hest))/np.max(np.abs(H)) )




# #--------------------------------------------
# HH = np.random.rand(N, N)
# HH = HH.T.dot(HH) + 10*np.eye(N)

# print(min(np.linalg.eigvals(HH)))

# f2 = f2_fact(HH)
# Hest_2 = Hess_est(x, h, f2)
# print("\n\n----------------------HH est")
# print(Hest_2 - HH)


# # print("HH")
# # print(HH)

# ipdb.set_trace()
