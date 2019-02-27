"""python example of using libdl1qc_dct.so with ctypes.
This is a lot easier than writing a custom wrapper.

https://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""
import matplotlib.pyplot as plt
import json
import time
import numpy as np
from ctypes import Structure
from ctypes import c_int
from ctypes import c_double
from ctypes import POINTER, byref


class LBResult(Structure):
    _fields_ = [("l1", c_double),
                ("total_newton_iter", c_int),
                ("total_cg_iter", c_int),
                ("status", c_int)]

    def __init__(self):
        self.l1 = 0
        self.total_newton_iter = 0
        self.total_cg_iter = 0
        self.status = 0

    def __repr__(self):
        return ""                                                 \
            "l1                : %f\n" % self.l1 +                \
            "total_newton_iter : %d\n" % self.total_newton_iter + \
            "total_cg_iter     : %d\n" % self.total_cg_iter +     \
            "status            : %d\n" % self.status


class l1qc_dct_params(Structure):
    _fields_ = [("epsilon", c_double),
                ("mu", c_double),
                ("lbtol", c_double),
                ("tau", c_double),
                ("lbiter", c_int),
                ("newton_tol", c_double),
                ("newton_max_iter", c_int),
                ("verbose", c_int),
                ("l1_tol", c_double),
                ("cgtol", c_double),
                ("cgmaxiter", c_int),
                ("warm_start_cg", c_int)]

    def __init__(self, epsilon=0.1, mu=10, lbtol=1e-3,
                 tau=0, lbiter=0,
                 newton_tol=1e-3, newton_max_iter=50,
                 verbose=2, l1_tol=0, cgtol=1e-8,
                 cgmaxiter=200, warm_start_cg=0):

        self.epsilon = epsilon
        self.mu = mu
        self.lbtol = lbtol
        self.tau = tau
        self.lbiter = lbiter
        self.newton_tol = newton_tol
        self.newton_max_iter = newton_max_iter
        self.verbose = verbose
        self.l1_tol = l1_tol
        self.cgtol = cgtol
        self.cgmaxiter = cgmaxiter
        self.warm_start_cg = warm_start_cg

    def __repr__(self):
        return "" + \
            "epsilon         : %f\n" % self.epsilon +         \
            "mu              : %f\n" % self.mu +              \
            "lbtol           : %f\n" % self.lbtol +           \
            "tau             : %f\n" % self.tau +             \
            "lbiter             : %f\n" % self.lbiter +       \
            "newton_tol      : %f\n" % self.newton_tol +      \
            "newton_max_iter : %d\n" % self.newton_max_iter + \
            "verbose         : %d\n" % self.verbose +         \
            "l1_tol          : %f\n" % self.l1_tol +          \
            "cgtol           : %f\n" % self.cgtol +           \
            "cgmaxiter       : %d\n" % self.cgmaxiter +       \
            "warm_start_cg   : %d" % self.warm_start_cg


def l1qc_dct(eta_0, b, pix_idx, opts, lib_dir=""):
    import numpy.ctypeslib as npct

    array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
    array_1d_int = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')

    libl1c = npct.load_library(lib_dir+"/libl1qc_dct.so", ".")

    libl1c.l1qc_dct.restype = np.int32
    libl1c.l1qc_dct.argtypes = [c_int, array_1d_double,
                                c_int, array_1d_double, array_1d_int,
                                l1qc_dct_params, POINTER(LBResult)]

    N = len(eta_0)
    M = len(b)
    eta_out = np.zeros(N)
    lb_res = LBResult()
    libl1c.l1qc_dct(N, eta_out, M, b, pix_idx, opts, byref(lb_res))

    return (eta_out, lb_res)


def remove_ticks(ax):
    ax.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',       # both major and minor ticks are affected
        bottom=False,       # ticks along the bottom edge are off
        top=False,          # ticks along the top edge are off
        left=False,         # ticks along left edge are off
        labelbottom=False,  # labels along the bottom edge are off
        labelleft=False,    # labels along the left edge are off
        labelright=False)   # labels along the right edge are off


def dct_example(verbose=2, fpath='example_img_data.json', plot=False, lib_dir=""):
    with open(fpath) as json_data:
        d = json.load(json_data)

        # pix_idx is a index set where we sampled. x_orig is a vector  of the
        # original image.
        pix_idx = np.array(d['pix_idx'], dtype=np.int32, ndmin=1)
        x_orig = np.array(d['x_orig'], ndmin=1)

        x_masked = np.zeros(len(x_orig))
        x_masked[pix_idx] = x_orig[pix_idx]
        b = x_orig[pix_idx]

        # Options for l1qc_dct()
        opts = l1qc_dct_params(epsilon=0.1,   # eps
                               mu=10,    # mu
                               lbtol=1e-3,  # lbtol
                               newton_tol=1e-3,  # newton_tol
                               newton_max_iter=50,    # newton_max_iter
                               verbose=verbose,     # verbose
                               l1_tol=1e-5,  # l1_tol
                               cgtol=1e-8,  # cgtol
                               cgmaxiter=200,   # cgmaxiter
                               warm_start_cg=0)     # warm_start_cg

        # Call the library wrapper.
        time0 = time.time()
        x_recon, lb_result = l1qc_dct(x_orig, b, pix_idx, opts, lib_dir=lib_dir)
        time_total = time.time() - time0

        print("Total python time: %f" % time_total)
        # Turn the vectors back into matrices so we can show them as an image.
        N = int(np.sqrt(len(x_orig)))

        X_orig_mat = np.reshape(x_orig, (N, N))
        X_recon_mat = np.reshape(x_recon, (N, N))
        X_masked_mat = np.reshape(x_masked, (N, N))
        if plot:
            plt.figure(num=1, figsize=(12, 4))

            ax1 = plt.subplot(131)
            remove_ticks(ax1)
            ax1.set_title("Original Image (CS-20ng grating)")
            ax1.imshow(X_orig_mat, cmap='gray')
            remove_ticks(ax1)

            ax2 = plt.subplot(132)
            ax2.imshow(X_masked_mat, cmap='gray')
            ax2.set_title("Subsampled image")
            remove_ticks(ax2)

            ax3 = plt.subplot(133)
            ax3.set_title("Reconstruction")
            ax3.imshow(X_recon_mat, cmap='gray')
            remove_ticks(ax3)

            plt.show()


if __name__ == '__main__':
    import sys
    import os
    plot = False

    if len(sys.argv) >= 2:
        fpath = sys.argv[1]
    else:
        fpath = "example_img_data.json"

    if len(sys.argv) >= 3:
        verbose = int(sys.argv[2])
    else:
        verbose = 1

    if len(sys.argv) >= 4:
        plot = True

    lib_dir = os.getenv("LIB_DIR")
    if lib_dir is None:
        lib_dir = ""

    dct_example(verbose=verbose, fpath=fpath, plot=plot, lib_dir=lib_dir)
