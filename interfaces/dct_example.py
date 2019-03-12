"""python example of using libdl1qc_dct.so with ctypes.
This is a lot easier than writing a custom wrapper.

https://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html
"""
import matplotlib.pyplot as plt
import json
import time
import numpy as np


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


def dct_example(verbose=2, fpath='example_img_data.json', plot=False):
    import l1cPy

    with open(fpath) as json_data:
        d = json.load(json_data)

    # pix_idx is a index set where we sampled. x_orig is a vector  of the
    # original image.
    pix_idx = np.array(d['pix_idx'], dtype=np.int32, ndmin=1)
    x_orig = np.array(d['x_orig'], ndmin=1)

    x_masked = np.zeros(len(x_orig))
    x_masked[pix_idx] = x_orig[pix_idx]
    b = x_orig[pix_idx]
    N = int(np.sqrt(len(x_orig)))

    # Call the library wrapper.
    time0 = time.time()
    x_recon, p = l1cPy.l1qc_dct(N*N, 1, b, pix_idx, epsilon=0.01,
                                   l1_tol=0, cgmaxiter=200, verbose=verbose)

    # x_recon, lb_result = l1qc_dct(x_orig, b, pix_idx, opts, lib_dir=lib_dir)
    time_total = time.time() - time0

    print("Total python time: %f" % time_total)
    # Turn the vectors back into matrices so we can show them as an image.

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
    import argparse

    description = 'Run the dct CS example solved via l1qc.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('fpath', nargs='?', action='store',
                        default="example_img_data.json")
    parser.add_argument('--plot', dest='plot', default=False, type=bool)
    parser.add_argument('--verbose', dest='verbose', default=0, type=int)

    args = parser.parse_args()

    plot = args.plot
    fpath = args.fpath
    verbose = args.verbose

    interface_dir = os.getenv("L1C_INTERFACE_DIR")
    src_lib_dir = os.getenv("L1C_SRC_LIB_DIR")

    if interface_dir is not None:
        sys.path.append(interface_dir)

    if src_lib_dir is not None:
        os.environ['PATH'] = src_lib_dir+";" + os.environ['PATH']
        lib_dir = src_lib_dir

    dct_example(verbose=verbose, fpath=fpath, plot=plot)
