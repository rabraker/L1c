import json
import L1c
import numpy as np
import matplotlib.pyplot as plt


def TV_denoise(verbose=2, fpath='example_img_data.json', plot=False, lib_dir="."):
    import dct_example

    with open(fpath) as json_data:
        d = json.load(json_data)

    x_orig = np.array(d['x_orig'], ndmin=1)
    n = int(np.sqrt(len(x_orig)))
    m = n
    # ipdb.set_trace()
    np.random.seed(0)
    x_noisy = x_orig + np.random.rand(n*m)

    # Call the library wrapper.
    mu = 5
    maxiter = 100
    max_jac_iter = 1
    x_recon = L1c.BRAnsTV(x_noisy, n, m, mu, maxiter,
                          max_jac_iter, lib_dir=lib_dir)

    # Turn the vectors back into matrices so we can show them as an image.

    X_orig_mat = np.reshape(x_orig, (n, m))
    X_noisy_mat = np.reshape(x_noisy, (n, m))
    X_recon_mat = np.reshape(x_recon, (n, m))
    if plot:
        plt.figure(num=1, figsize=(12, 4))

        ax1 = plt.subplot(131)
        ax1.set_title("Original Image (CS-20ng grating)")
        ax1.imshow(X_orig_mat, cmap='gray')
        dct_example.remove_ticks(ax1)

        ax2 = plt.subplot(132)
        ax2.imshow(X_noisy_mat, cmap='gray')
        ax2.set_title("Noisy image")
        dct_example.remove_ticks(ax2)

        ax3 = plt.subplot(133)
        ax3.set_title("Reconstruction")
        ax3.imshow(X_recon_mat, cmap='gray')
        dct_example.remove_ticks(ax3)

        plt.show()


fpath = "../test/test_data/example_img_data.json"
libdir = "../build-fftw-ob/src/.libs"
TV_denoise(verbose=0, fpath=fpath, plot=True, lib_dir=libdir)
