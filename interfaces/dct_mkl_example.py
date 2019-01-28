"""python example of using libdl1qc_dct.so with ctypes.
This is a lot easier than writing a custom wrapper.

https://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""
import matplotlib.pyplot as plt
import json
import numpy as np
from ctypes import Structure
from ctypes import c_int
from ctypes import c_double

class l1qc_dct_params(Structure):
    _fields_ = [("epsilon", c_double),
                ("mu", c_double),
                ("lbtol", c_double),
                ("newton_tol", c_double),
                ("newton_max_iter", c_int),
                ("verbose", c_int),
                ("l1_tol", c_double),
                ("cgtol", c_double),
                ("cgmaxiter", c_int),
                ("warm_start_cg", c_int)]


def l1qc_dct(eta_0, b, pix_idx, opts):
    import numpy.ctypeslib as npct

    array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
    array_1d_int = npct.ndpointer(dtype=np.int64, ndim=1, flags='CONTIGUOUS')

    libl1c = npct.load_library("libl1qc_dct", ".")

    libl1c.l1qc_dct.restype = np.int32
    libl1c.l1qc_dct.argtypes = [c_int, array_1d_double,
                                c_int, array_1d_double, array_1d_int,
                                l1qc_dct_params]

    N = len(eta_0)
    M = len(b)
    eta_out = np.zeros(N)
    libl1c.l1qc_dct(N, eta_out, M, b, pix_idx, opts)

    return eta_out


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


def dct_mkl_example():
    with open('example_img_data.json') as json_data:
        d = json.load(json_data)

        # pix_idx is a index set where we sampled. x_orig is a vector  of the
        # original image.
        pix_idx = np.array(d['pix_idx'], dtype=np.int64, ndmin=1)
        x_orig = np.array(d['x_orig'], ndmin=1)

        x_masked = np.zeros(len(x_orig))
        x_masked[pix_idx] = x_orig[pix_idx]
        b = x_orig[pix_idx]

        # Options for l1qc_dct()
        opts = l1qc_dct_params(0.1,   # eps
                               10,    # mu
                               1e-3,  # lbtol
                               1e-3,  # newton_tol
                               50,    # newton_max_iter
                               2,     # verbose
                               1e-5,  # l1_tol
                               1e-8,  # cgtol
                               200,   # cgmaxiter
                               0)     # warm_start_cg
        # Call the library wrapper.
        x_recon = l1qc_dct(x_orig, b, pix_idx, opts)

        # Turn the vectors back into matrices so we can show them as an image.
        N = int(np.sqrt(len(x_orig)))

        X_orig_mat = np.reshape(x_orig, (N, N))
        X_recon_mat = np.reshape(x_recon, (N, N))
        X_masked_mat = np.reshape(x_masked, (N, N))

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
    dct_mkl_example()
