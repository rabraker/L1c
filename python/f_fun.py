import numpy as np
from scipy.fftpack import dct
import ipdb
import json
import codecs


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


def build_h11p_data(A, pix_idx):
    h11p_path = test_data_path+"/hp11_fun_data.json"

    N = A.shape[1]
    M = A.shape[0]

    sigx = np.random.randn(N, 1)
    z = np.random.randn(N, 1)
    r = np.random.randn(M, 1)

    atr = A.T.dot(r)
    fe = np.random.rand(1)[0]
    fe = -np.abs(fe)

    h11p_z = H11p_fun(sigx, fe, A, atr, z)

    data = {'z': z,
            'atr': atr,
            'fe': fe,
            'sigx': sigx,
            'y_exp': h11p_z,
            'pix_idx': pix_idx}

    data = jsonify(data)
    save_json(data, h11p_path)


def build_feval_data(x, u, b, Amat, epsilon, tau):
    feval_path = test_data_path+"/f_eval_data.json"

    xu = np.vstack((x, u))
    r, f, fe, fu1, fu2 = f_fun(xu, b, Amat, epsilon, tau)

    data = {'x': x,
            'u': u,
            'r': r,
            'tau': tau,
            'epsilon': epsilon,
            'fu1_exp': fu1,
            'fu2_exp': fu2,
            'fe_exp': fe,
            'f_exp': f}

    data = jsonify(data)

    save_json(data, feval_path)


def find_max_step(dx, du, Adx, fu1, fu2, r, epsilon):
    # ipdb.set_trace()
    idx_fu1, = np.where((dx-du)[:, 0] > 0)
    idx_fu2, = np.where((-dx-du)[:, 0] > 0)

    aqe = Adx.T.dot(Adx)[0][0]
    bqe = 2*r.T.dot(Adx)[0][0]
    cqe = (r.T.dot(r) - epsilon**2)[0][0]

    smax_f1 = np.min(-fu1[idx_fu1] / (dx[idx_fu1] - du[idx_fu1]))
    smax_f2 = np.min(-fu2[idx_fu2] / (-dx[idx_fu2] - du[idx_fu2]))
    smax_quad = (-bqe + np.sqrt(bqe**2 - 4*aqe*cqe)) / (2*aqe)

    print(smax_f1)
    print(smax_f2)
    print(smax_quad)
    smax = np.min([1, smax_f1, smax_f2, smax_quad])

    return smax * 0.99


def build_find_smax_data(dx, du, Adx, fu1, fu2, r, epsilon):
    smax_path = test_data_path+"/find_max_step_data.json"

    print(smax_path)
    smax = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon)

    data = {'dx': dx,
            'du': du,
            'Adx': Adx,
            'fu1': fu1,
            'fu2': fu2,
            'r': r,
            'epsilon': epsilon,
            'smax': smax}

    data = jsonify(data)
    save_json(data, smax_path)


def newton_init(x0, lbtol, mu):
    N = x0.shape[0]
    u = (0.95)*np.abs(x0) + (0.10)*np.max(np.abs(x0))
    ipdb.set_trace()
    tmp = (2.0*N+1.0)/np.sum(np.abs(x0))
    tau = np.max((tmp, 1))
    lbiter = np.ceil((np.log(2.0*N+1.0)-np.log(lbtol)-np.log(tau))/np.log(mu))

    return u, tau, lbiter


def build_newton_init_data(x0, lbtol, mu, epsilon, b, pix_idx):
    newton_init_path = test_data_path+"/newton_init_data.json"

    u, tau, lbiter = newton_init(x0, lbtol, mu)

    data = {'x': x,
            'u': u,
            'b': b,
            'lbtol': lbtol,
            'mu': mu,
            'tau': tau,
            'epsilon': epsilon,
            'lbiter': lbiter,
            'pix_idx': pix_idx}

    data = jsonify(data)
    save_json(data, newton_init_path)


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

    ntgz = 1/fu1 - 1/fu2 + (1/fe) * atr
    ntgu = -tau - 1/fu1 - 1/fu2
    gradf = -(1/tau)*np.hstack((ntgz, ntgu))

    return gradf


def H11p_fun(sigx, fe, A, atr, z):
    # Need the slash here, python will not continue the line on its own.
    h11p_z = sigx*z + (-1.0/fe) * A.T.dot(A.dot(z)) + \
             (1/fe**2) * (atr.T.dot(z))*atr

    return h11p_z




np.random.seed(0)
N = 64
tau = 10
mu = 10
lbtol = 1e-3
epsilon = .1
Amat, Atmat, pix_idx = build_Amats(N)

M = Amat.shape[0]
xx = np.random.randn(N, 1)
b = Amat.dot(xx) + np.random.randn(M, 1) * 0.015
x = Atmat.dot(b) + np.random.randn(N, 1) * 0.0001
u = (0.99)*np.abs(x) + (0.10)*np.max(np.abs(x))

dx = np.random.randn(N, 1)
du = np.random.randn(N, 1)
Adx = Amat.dot(dx)

r, f, fe, fu1, fu2 = f_fun(np.vstack((x, u)), b, Amat, epsilon, tau)
smax = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon)

build_find_smax_data(dx, du, Adx, fu1, fu2, r, epsilon)

build_newton_init_data(x, lbtol, mu, epsilon, b, pix_idx)


# build_h11p_data(Amat, pix_idx)
# build_feval_data(x, u, b, Amat, epsilon, tau)

print("||Ax-b||")
print(np.linalg.norm(Amat.dot(x) - b) )
ipdb.set_trace()


# xu = np.hstack((x, u))
# r, f, fe, fu1, fu2 = f_fun(xu, b, Amat, epsilon, tau)

# f_fun_handle = f_fun_fact(b, Amat, epsilon, tau)



# gf_true = gradf(Amat, r, fu1, fu2, fe, tau)




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
