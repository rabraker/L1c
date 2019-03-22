#!/usr/bin/env python3
import os
import numpy as np
import L1cTestDataUtils as TDU
from numpy.random import seed


def build_matrix_xfms_test_data(test_data_path, N, M):
    # Build and save test data for matrix_transforms.c
    #

    seed(2)
    E_idx = np.random.randint(0, N, M)
    A = np.random.rand(N, N)[E_idx, :]

    x = np.random.rand(N)*10
    y = np.random.rand(M)*10

    Ax = A.dot(x)
    Aty = A.T.dot(y)

    AtAx = A.T.dot(A.dot(x))

    data = {'Nrow': N,
            'Mcol': M,
            'A': A,
            'x': x,
            'y': y,
            'Ax': Ax,
            'Aty': Aty,
            'AtAx': AtAx}

    data = TDU.jsonify(data)
    TDU.save_json(data, test_data_path)


if __name__ == "__main__":
    srcdir = os.getenv("srcdir")
    if srcdir is None:
        srcdir = "."

    data_dir = srcdir+"/test_data"

    N = 50
    M = 20

    small_matrix_path = data_dir+"/matrix_xfm_small_data.json"
    build_matrix_xfms_test_data(small_matrix_path, N, M)
