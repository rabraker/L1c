import json
import time
import numpy as np

import l1c_pyplot_utils as lpu
import l1c_py_init_path as lpip


def breg_anisTV_example(fpath='example_img_data.json', plot=False):
    import _l1cPy_module as l1cPy

    with open(fpath) as json_data:
        d = json.load(json_data)

    x_orig = np.array(d['x_orig'], ndmin=1)
    n = int(np.sqrt(len(x_orig)))
    m = n

    np.random.seed(0)
    x_noisy = x_orig + np.random.rand(n*m)
    X_noisy_mat = np.reshape(x_noisy, (n, m))

    start = time.perf_counter()
    Xclean_mat = l1cPy.breg_anistropic_TV(X_noisy_mat,
                                          max_iter=100, max_jac_iter=1,
                                          tol=0.001, mu=5)
    end = time.perf_counter()
    print("(breg_anistropic_TV) Python time = %f" % (end - start))

    if plot:
        import matplotlib.pyplot as plt

        plt.figure(num=1, figsize=(8, 4))

        ax1 = plt.subplot(121)
        lpu.remove_ticks(ax1)
        ax1.set_title("Original Image")
        ax1.imshow(X_noisy_mat, cmap='gray')
        lpu.remove_ticks(ax1)

        ax2 = plt.subplot(122)
        ax2.imshow(Xclean_mat, cmap='gray')
        ax2.set_title("Anistropic TV denoised")
        lpu.remove_ticks(ax2)

        plt.show()


if __name__ == "__main__":
    fpath = "example_img_data.json"
    lpip.add_lib_path()
    breg_anisTV_example(fpath=fpath, plot=True)
