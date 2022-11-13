#!/usr/bin/env python3
"""python example of using nesta_dctTv
"""

import json
import time

import l1c_py_init_path as lpip
import l1c_pyplot_utils as lpu
import numpy as np


def dct_example(verbose=2, fpath="example_img_data.json", plot=False):
    import _l1cPy_module as l1cPy

    with open(fpath) as json_data:
        d = json.load(json_data)

    # pix_idx is a index set where we sampled. x_orig is a vector  of the
    # original image.
    pix_idx = np.array(d["pix_idx"], dtype=np.int32, ndmin=1)
    x_orig = np.array(d["x_orig"], ndmin=1)

    x_masked = np.zeros(len(x_orig))
    x_masked[pix_idx] = x_orig[pix_idx]
    b = x_orig[pix_idx]
    N = int(np.sqrt(len(x_orig)))

    alpha_v = 0.5
    alpha_h = 0.1
    tol = 1e-5
    mu = 1e-5
    time0 = time.perf_counter()
    x_recon, status = l1cPy.nesta_dctTV(
        N, N, b, pix_idx, alpha_v=alpha_v, alpha_h=alpha_h, mu=mu, tol=tol
    )

    time_total = time.perf_counter() - time0

    print("(l1qc_dct) Total python time: %f" % time_total)

    if plot:
        import matplotlib.pyplot as plt

        # Turn the vectors back into matrices so we can show them as an image.
        X_orig_mat = np.reshape(x_orig, (N, N))
        X_recon_mat = np.reshape(x_recon, (N, N))
        X_masked_mat = np.reshape(x_masked, (N, N))

        plt.figure(num=1, figsize=(12, 4))

        ax1 = plt.subplot(131)
        lpu.remove_ticks(ax1)
        ax1.set_title("Original Image (CS-20ng grating)")
        ax1.imshow(X_orig_mat, cmap="gray")
        lpu.remove_ticks(ax1)

        ax2 = plt.subplot(132)
        ax2.imshow(X_masked_mat, cmap="gray")
        ax2.set_title("Subsampled image")
        lpu.remove_ticks(ax2)

        ax3 = plt.subplot(133)
        ax3.set_title("Reconstruction")
        ax3.imshow(X_recon_mat, cmap="gray")
        lpu.remove_ticks(ax3)

        plt.show()


if __name__ == "__main__":
    import argparse
    import os
    import sys

    description = "Run the dct CS example solved via l1qc."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "fpath", nargs="?", action="store", default="example_img_data.json"
    )
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--verbose", dest="verbose", default=0, type=int)

    args = parser.parse_args()

    plot = args.plot
    fpath = args.fpath
    verbose = args.verbose

    print("plot = %s" % plot)
    print("verbose = %d" % verbose)
    lpip.add_lib_path()

    dct_example(verbose=verbose, fpath=fpath, plot=plot)
