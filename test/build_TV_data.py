#!/usr/bin/env python3
import L1cTestDataUtils as TDU
import numpy as np
import numpy.testing as npt


def Dy_as_diff(x_vec, n, m):
    zro_row = np.zeros((1, M))
    X = np.reshape(x_vec, (n, m), order="C")

    DyX = np.vstack((np.diff(X, axis=0), zro_row))
    DyX = np.reshape(DyX, (n * m, 1), order="C")
    return DyX


def Dx_as_diff(x_vec, n, m):
    zro_col = np.zeros((N, 1))
    X = np.reshape(x_vec, (n, m), order="C")

    DxX = np.hstack((np.diff(X, axis=1), zro_col))
    DxX = np.reshape(DxX, (n * m, 1), order="C")
    return DxX


def DxMatRep(n, m):
    Dx_kernel = -np.eye(m) + np.hstack((np.zeros((m, 1)), np.eye(m, m - 1)))
    Dx_kernel[-1, -1] = 0
    In = np.eye(n)
    Dx = np.kron(In, Dx_kernel)

    return Dx


def DyMatRep(n, m):
    Dy_kernel = -np.eye(n) + np.hstack((np.zeros((n, 1)), np.eye(n, n - 1)))
    Dy_kernel[-1, -1] = 0
    Im = np.eye(m)
    Dy = np.kron(Dy_kernel, Im)

    return Dy


def Dy(x_vec, n, m):
    """
    [[-1. -0. -0.  1.  0.  0.  0.  0.  0.  0.  0.  0.]
     [-0. -1. -0.  0.  1.  0.  0.  0.  0.  0.  0.  0.]
     [-0. -0. -1.  0.  0.  1.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0. -1. -0. -0.  1.  0.  0.  0.  0.  0.]
     [ 0.  0.  0. -0. -1. -0.  0.  1.  0.  0.  0.  0.]
     [ 0.  0.  0. -0. -0. -1.  0.  0.  1.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0. -1. -0. -0.  1.  0.  0.]
     [ 0.  0.  0.  0.  0.  0. -0. -1. -0.  0.  1.  0.]
     [ 0.  0.  0.  0.  0.  0. -0. -0. -1.  0.  0.  1.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]
    """

    dy = np.zeros((n * m, 1)) + 1  # Ensure everything is set right.

    for row in range(0, n - 1):  # n-2 inclusive
        for col in range(row * m, (row + 1) * m):  # (row+1)*m-1 inclusive
            i = col
            dy[i] = x_vec[i + m] - x_vec[i]

    for i in range((n - 1) * m, n * m):
        dy[i] = 0
    return dy


def Dx(x_vec, n, m):
    """
        4 by 3
    [[-1.  1.  0. -0.  0.  0. -0.  0.  0. -0.  0.  0.]
     [ 0. -1.  1.  0. -0.  0.  0. -0.  0.  0. -0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [-0.  0.  0. -1.  1.  0. -0.  0.  0. -0.  0.  0.]
     [ 0. -0.  0.  0. -1.  1.  0. -0.  0.  0. -0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [-0.  0.  0. -0.  0.  0. -1.  1.  0. -0.  0.  0.]
     [ 0. -0.  0.  0. -0.  0.  0. -1.  1.  0. -0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [-0.  0.  0. -0.  0.  0. -0.  0.  0. -1.  1.  0.]
     [ 0. -0.  0.  0. -0.  0.  0. -0.  0.  0. -1.  1.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]]
    """

    dx = np.zeros((n * m, 1)) + 1  # Ensure everything is set right.

    for row in range(0, n):  # n-1 inclusive
        for col in range(row * m, (row + 1) * m - 1):  # (row+1)*m - 1 inclusive
            i = col  # row*m +col + 1
            dx[i] = x_vec[i + 1] - x_vec[i]

        dx[i + 1] = 0
    return dx


def DyTDy(A_vec, alpha, n, m):
    """
     Given an m by n matrix A, computes lambda*(Del_y^T*Del_y)*A.
     We assume A is stored in the 1-D vector A, in row major order.
     For a 3 x 4 matrix, Del_y^T*Del_y has the matrix representation:

     1     0     0    -1     0     0     0     0     0     0     0     0
     0     1     0     0    -1     0     0     0     0     0     0     0
     0     0     1     0     0    -1     0     0     0     0     0     0
    -1     0     0     2     0     0    -1     0     0     0     0     0
     0    -1     0     0     2     0     0    -1     0     0     0     0
     0     0    -1     0     0     2     0     0    -1     0     0     0
     0     0     0    -1     0     0     2     0     0    -1     0     0
     0     0     0     0    -1     0     0     2     0     0    -1     0
     0     0     0     0     0    -1     0     0     2     0     0    -1
     0     0     0     0     0     0    -1     0     0     1     0     0
     0     0     0     0     0     0     0    -1     0     0     1     0
     0     0     0     0     0     0     0     0    -1     0     0     1
    """
    dytdy = A_vec * 0 + 1  # Ensure everything is set right.

    # First, do the diagonal and upper diagonal
    for i in range(0, m * n):
        if i < m:
            Ai_min_m = 0
        else:
            Ai_min_m = A_vec[i - m]

        if (i < m) or (i > m * n - m - 1):
            D_ii_Ai = A_vec[i]
        else:
            D_ii_Ai = 2 * A_vec[i]

        if i < n * m - m:
            Ai_p_m = A_vec[i + m]
        else:
            Ai_p_m = 0

        dytdy[i] = (alpha**2) * (-Ai_min_m + D_ii_Ai - Ai_p_m)

    return dytdy


def DxTDx(A_vec, alpha, n, m):
    """
    [[ 1. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [-1.  2. -1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0. -1.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  1. -1.  0.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0. -1.  2. -1.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0. -1.  1.  0.  0.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  1. -1.  0.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0. -1.  2. -1.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0. -1.  1.  0.  0.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  1. -1.  0.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  2. -1.]
     [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0. -1.  1.]]
    """

    dxtdx = A_vec * 0 + 1  # Ensure everything is set right.

    for row in range(0, n):
        for col in range(0, m):
            i = row * m + col

            if col == 0 or col == m - 1:
                D_ii = A_vec[i]
            else:
                D_ii = 2.0 * A_vec[i]

            if col == 0:
                Ai_m1 = 0
            else:
                Ai_m1 = A_vec[i - 1]

            if col == m - 1 or i == n * m - 1:
                Ai_p1 = 0
            else:
                Ai_p1 = A_vec[i + 1]

            dxtdx[i] = (alpha**2) * (-Ai_m1 + D_ii - Ai_p1)

    return dxtdx


def check_Dx(A, N, M):
    a_vec = np.reshape(A, (N * M, 1), order="C")
    Dx_vec_exp = Dx_as_diff(a_vec, N, M)
    DxMat = DxMatRep(N, M)
    Dx_loop = Dx(a_vec, N, M)

    npt.assert_array_almost_equal(Dx_loop, Dx_vec_exp, decimal=14)
    npt.assert_array_almost_equal(Dx_vec_exp, DxMat.dot(a_vec), decimal=14)
    npt.assert_array_almost_equal(Dx_loop, DxMat.dot(a_vec), decimal=14)


def check_Dy(A, N, M):
    a_vec = np.reshape(A, (N * M, 1), order="C")
    Dy_vec_exp = Dy_as_diff(a_vec, N, M)
    Dy_loop = Dy(a_vec, N, M)
    DyMat = DyMatRep(N, M)

    npt.assert_array_almost_equal(Dy_vec_exp, Dy_loop, decimal=14)

    npt.assert_array_almost_equal(Dy_vec_exp, DyMat.dot(a_vec), decimal=14)
    npt.assert_array_almost_equal(Dy_loop, DyMat.dot(a_vec), decimal=14)


def check_DxTDx(A, N, M):
    alpha = 2.5436
    a_vec = np.reshape(A, (N * M, 1), order="C")
    dxTdx = DxTDx(a_vec, alpha, N, M)

    DxMat = DxMatRep(N, M)

    dxTdx_exp = alpha * DxMat.T.dot(alpha * DxMat.dot(a_vec))
    npt.assert_array_almost_equal(dxTdx, dxTdx_exp, decimal=12)


def check_DyTDy(A, N, M):
    alpha = 2.5436
    a_vec = np.reshape(A, (N * M, 1), order="C")
    dyTdy = DyTDy(a_vec, alpha, N, M)

    DyMat = DyMatRep(N, M)
    dyTdy_exp = alpha * DyMat.T.dot(alpha * DyMat.dot(a_vec))
    npt.assert_array_almost_equal(dyTdy, dyTdy_exp, decimal=12)


# To build the laplacian, we need a matrix representation of Dy^T and Dx^T.
# It is easy to take the gradient in y or x direction with a diff along the
# two axes. It is not so obvious what the adjoint of the diff operation is.
# The strategy here is to build
#   1. the niave diff in x and y (Dx_as_diff, Dy_as_Diff)
#   2. Build a matrix representiation of Dx and Dy.
#   3. Assert that the give the same result.
#   4. Assert that Dx*x_vec etc give the same result as a for-loop based
#      implementation (which is similar to the c-code).


def build_TV_data(N, M):
    alpha = 2.67
    # A = np.random.rand(N, M)
    A = np.random.randint(1, high=100, size=(N, M))
    A = np.double(A)
    check_Dx(A, N, M)
    check_Dy(A, N, M)
    check_DxTDx(A, N, M)
    check_DyTDy(A, N, M)

    a_vec = np.reshape(A, (N * M, 1), order="C")
    DyMat = DyMatRep(N, M)
    DxMat = DxMatRep(N, M)

    dxA = alpha * DxMat.dot(a_vec)
    dyA = alpha * DyMat.dot(a_vec)

    dxTA = alpha * DxMat.T.dot(a_vec)
    dyTA = alpha * DyMat.T.dot(a_vec)

    dxTdxA = alpha * DxMat.T.dot(alpha * DxMat.dot(a_vec))
    dyTdyA = alpha * DyMat.T.dot(alpha * DyMat.dot(a_vec))
    data = {
        "A": a_vec,
        "DxA": dxA,
        "DyA": dyA,
        "DxTA": dxTA,
        "DyTA": dyTA,
        "DxTDxA": dxTdxA,
        "DyTDyA": dyTdyA,
        "alpha": alpha,
        "N": N,
        "M": M,
        "NM": N * M,
    }

    data = TDU.jsonify(data)

    return data


if __name__ == "__main__":
    data_dir = TDU.data_dir()

    N = 16
    M = 16
    data = build_TV_data(N, M)
    fname_square = data_dir + "/TV_data_square_small.json"
    TDU.save_json(data, fname_square)

    N = 16
    M = 13
    data = build_TV_data(N, M)
    fname_square = data_dir + "/TV_data_tall_skinny.json"
    TDU.save_json(data, fname_square)

    N = 13
    M = 16
    data = build_TV_data(N, M)
    fname_square = data_dir + "/TV_data_short_wide.json"
    TDU.save_json(data, fname_square)

# """
# typedef struct ImDiffData{
#   double *A;
#   double *DxA_exp;
#   double *DyA_exp;

#   double *DxTA_exp;
#   double *DyTA_exp;

#   double *DxTDxA_exp;
#   double *DyTDyA_exp;
#   :
#   :
#   :
#   int N;
#   int M;
#   int NM;
#   char fpath[256];
# }ImgDiffData;
# """
