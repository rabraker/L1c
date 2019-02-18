import numpy as np
import json
import codecs

from numpy.random import seed, rand
from scipy.fftpack import dct


def save_json(data, fname):
    json.dump(data, codecs.open(fname, 'w', encoding='utf-8'),
              separators=(',', ':'), sort_keys=True, indent=4)


def Adct_factory(pix_idx):
    def Adct(x):
        y = dct(x, type=3, norm='ortho')
        y = y[pix_idx]
        return y
    return Adct


def Atdct_factory(pix_idx, N):
    def Atdct(y):
        x = np.zeros(N)
        x[pix_idx] = y
        x = dct(x, type=2, norm='ortho')
        return x

    return Atdct


def build_dct_rand_test_data(fname, pix_idx, Nx):

    seed(0)

    eta_vec = rand(Nx)

    Adct = Adct_factory(pix_idx)
    Atdct = Atdct_factory(pix_idx, Nx)

    # Another small example with randomly generated data.

    # y = E*M*x
    y_vec = dct(eta_vec, norm='ortho', type=3)[pix_idx]

    EMx = Adct(eta_vec)
    MtEty = Atdct(y_vec)
    MtEt_EMx = Atdct(Adct(eta_vec))

    # --------------------- M^TE^T EM ----------------------
    # data = struct('x0', x(:)', 'x1', x_AtA(:)', 'pix_idx', pix_idx(:)'-1)
    data = {'x_in': eta_vec.tolist(),
            'y_in': y_vec.tolist(),
            'EMx': EMx.tolist(),
            'pix_idx': pix_idx.tolist(),
            'MtEty': MtEty.tolist(),
            'MtEt_EMx': MtEt_EMx.tolist()}

    save_json(data, fname)

    # savejson('', struct('x_in', img_vec(:)', 'y_in', y_vec(:)',...
    #   'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx', pix_idx(:)'-1), jopts)


Nx = 50
pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])

fname = 'dct_small.json'

# This checks that we get the right result for a pure dct, and in particular,
# that we got the scaling at the first element correct.
build_dct_rand_test_data(fname, pix_idx, Nx)
pix_idx = np.arange(50)
fname = 'dct_small_pure_dct.json'
build_dct_rand_test_data(fname, pix_idx, Nx)

# N = 50
# Ts = 1/50
# To = 0.25

# omega = 2*np.pi/To

# k = np.arrange(0, N)

# x = rand(N)

# pix_idx = np.array([2, 10, 15, 20, 25, 30, 35, 40, 45, 49])

# Adct = Adct_factory(pix_idx)
# Atdct = Atdct_factory(pix_idx, N)

# yy = Adct(x)

# data = {'x0': x.tolist(),
#         'x1': yy.tolist(),
#         'pix_idx': pix_idx.tolist()}

# save_json(data, EMx_path)

# # data = struct('x0', x(:)', 'x1', yy(:)', 'pix_idx', pix_idx(:)'-1)

# # --------------------- M^TE^T ----------------------
# x_AtA = Atdct(yy)


# # --------------------- M^TE^T EM ----------------------


# Now, Use data from an actual compressed sensing situation.

# img_dat = load(fullfile(data_root, 'test_image_data.mat'))
# xorig = img_dat.xorig
# pix_idx = img_dat.pix_idx
# N = length(xorig)

# A = @(x) L1qcTestData.Afun_dct(x, pix_idx)  # E*M
# At = @(x) L1qcTestData.Atfun_dct(x, pix_idx, N)  #E^T*M^T

# b = xorig(pix_idx)

# EMx = A(xorig)
# MtEty = At(b)

# MtEt_EMx = At(A(xorig))

# jopts.FileName = fullfile(data_root, 'dct_large.json')
# savejson('', struct('x_in', xorig(:)', 'y_in', b(:)',...
# 'EMx', EMx(:)', 'MtEty', MtEty(:)', 'MtEt_EMx', MtEt_EMx(:)', 'pix_idx',...
# pix_idx(:)'-1), jopts)
