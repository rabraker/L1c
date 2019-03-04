import numpy as np
import matplotlib.pyplot as plt


def grad_est(x0, h, fun):
    """
    Compute an estimate of the functional fun via forward finite
    difference.
    """
    N = len(x0)

    dx = np.zeros(N)
    ek = np.zeros(N)
    fx0 = fun(x0)
    for k in range(0, N):
        ek = ek*0
        ek[k] = 1
        dx[k] = (fun(x0 + ek*h) - fx0) / h

    return dx


def fun1(c):
    def fun(x):
        return np.dot(c, x)

    return fun


def l1c_fun(x, u, r, tau, epsilon):
    fu1 = x - u
    fu2 = -x - u
    fe = 0.5 * (r.dot(r) - epsilon**2)

    f = np.sum(u) - (1/tau) * (sum(np.log(-fu1)) +
                               sum(np.log(-fu2)) + np.log(-fe))

    return fu1, fu2, fe, f


np.random.rng(1)

N = 10  # problem size

x = np.random.randn(N)
u = (0.95)*np.abs(x) + 0.1*np.max(np.abs(x))

r = np.random.randn(N)*0.001


h = 0.1
x0 = np.random.rand(N)
c = np.ones(N)
# fun1 = x^2, d dfun1/dx = 2x

dx0_exp = c

fun = fun1(c)
ge = grad_est(x0, h, fun)

print(ge)
print(c)
