#include "config.h"

#include <math.h>

#include "cblas.h"
#include "l1c.h"

/*
  This file provides a set of functions for computing the discrete
  gradient and laplacian of an N by M image represented as a vector
  in \mathbb{R}^{NM}, which has been concatenated ROW WISE. N.B., that
  this is different from the standard vec(X) operation, which builds a vector
  column wise. However, representing the vector as a row-wise stacking
  conforms to the c convention better.

 */
/* Forward declarations */

void l1c_Dx(l1c_int n, l1c_int m, double alpha, double* X, double* Dx);

void l1c_DxT(l1c_int n, l1c_int m, double alpha, double* A, double* dxt);

void l1c_DxTDx(l1c_int n, l1c_int m, double alpha, double* A, double* dxtdx);

void l1c_DyT(l1c_int n, l1c_int m, double alpha, double* A, double* dyt);

void l1c_DyTDy(l1c_int n, l1c_int m, double alpha, double* A, double* dytdy);

void l1c_Dy(l1c_int n, l1c_int m, double alpha, double* X, double* Dy);

/*
  I've changed this up from the initial implementation
  (commit a98f800fa6144ad7e2c8bad78a047a95c738ecad and prior).
  This new scheme is faster, more readible I think, and doesnt
  depend on __builtin_assume_aligned, which is the main motivation.
  Because __builtin_assume_aligned makes it difficult/cumbersome
  to pull these operators into the overcomplete dictionary for nesta.
 */

/* Basic forward difference, with zero boundary condition at the
   far end

   [-1.  1   0   0  0]
   [ 0. -1.  1   0  0]
   [ 0.  0. -1.  1  0]
   [ 0.  0   0. -1  1]
   [ 0   0   0   0  0]
*/
static inline void ddiff_z(l1c_int n, double alpha, double* x, double* dx) {
  for (int i = 1; i < n; i++) {
    dx[i - 1] = alpha * (x[i] - x[i - 1]);
  }
  dx[n - 1] = 0;
}

/*
  Transpose of the operator defined by ddif_z
  [-1.  0   0   0  0]
  [ 1. -1.  0   0  0]
  [ 0.  1. -1.  0  0]
  [ 0.  0   1. -1  0]
  [ 0   0   0   1  0]
  dx[0] = -x[0]
  dx[1] = dx[0] - dx[1]
*/
static inline void ddiffT_z(l1c_int n, double alpha, double* x, double* dx) {
  dx[0] = -alpha * x[0];
  for (int i = 1; i < n; i++) {
    dx[i] = alpha * (x[i - 1] - x[i]);
  }
  dx[n - 1] = alpha * x[n - 2];
}

/*
  Second derivitive
  [[ 1. -1.  0.  0]
  [-1.  2. -1.  0]
  [ 0  -1   2  -1]
  [ 0.  0   -1  1]
*/
static inline void ddiff2_z(l1c_int n, double alpha, double* x, double* dx) {
  double alp2 = alpha * alpha;
  dx[0] = alp2 * (x[0] - x[1]);
  for (int i = 1; i < n - 1; i++) {
    dx[i] = alp2 * (-x[i - 1] + 2 * x[i] - x[i + 1]);
  }
  dx[n - 1] = alp2 * (x[n - 1] - x[n - 2]);
}

/*
  Implements the gradient in x direction. For a
  4 by 3 matrix X  concatenated row-wise, this is equivalent
  to multiplication by the following matrix:

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
 */
void l1c_Dx(
    l1c_int n, l1c_int m, double alpha, double* restrict X, double* restrict Dx) {
  /*
    We go down the rows of the matrix X, and compute the difference
    Dx = X[i,j+1] - X[i,j]
  */
  double *X_row = NULL, *Dx_row = NULL;
#pragma omp parallel for private(X_row, Dx_row)
  for (int i_row = 0; i_row < n; i_row++) {
    X_row = X + i_row * m;
    Dx_row = Dx + i_row * m;
    ddiff_z(m, alpha, X_row, Dx_row);
  }
}

/*
  4 by 3 case:
  [[-1.  0.  0. -0.  0.  0. -0.  0.  0. -0.  0.  0.]
  [ 1. -1.  0.  0. -0.  0.  0. -0.  0.  0. -0.  0.]
  [ 0.  1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.]
  [-0.  0.  0. -1.  0.  0. -0.  0.  0. -0.  0.  0.]
  [ 0. -0.  0.  1. -1.  0.  0. -0.  0.  0. -0.  0.]
  [ 0.  0.  0.  0.  1.  0.  0.  0.  0.  0.  0.  0.]
  [-0.  0.  0. -0.  0.  0. -1.  0.  0. -0.  0.  0.]
  [ 0. -0.  0.  0. -0.  0.  1. -1.  0.  0. -0.  0.]
  [ 0.  0.  0.  0.  0.  0.  0.  1.  0.  0.  0.  0.]
  [-0.  0.  0. -0.  0.  0. -0.  0.  0. -1.  0.  0.]
  [ 0. -0.  0.  0. -0.  0.  0. -0.  0.  1. -1.  0.]
  [ 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  1.  0.]]
 */
void l1c_DxT(
    l1c_int n, l1c_int m, double alpha, double* restrict X, double* restrict DxT) {
  double *X_row = NULL, *Dxt_row = NULL;
#pragma omp parallel for private(X_row, Dxt_row)
  for (int i_row = 0; i_row < n; i_row++) {
    X_row = X + i_row * m;
    Dxt_row = DxT + i_row * m;
    ddiffT_z(m, alpha, X_row, Dxt_row);
  }
}

/*
4 by 3:
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
 */
void l1c_DxTDx(
    l1c_int n, l1c_int m, double alpha, double* restrict X, double* restrict DxTDx) {

  double *X_row = NULL, *DxTDx_row = NULL;
#pragma omp parallel for private(X_row, DxTDx_row)
  for (int i_row = 0; i_row < n; i_row++) {
    X_row = X + i_row * m;
    DxTDx_row = DxTDx + i_row * m;
    ddiff2_z(m, alpha, X_row, DxTDx_row);
  }
}

/*
  ---------- PERFORMANCE NOTE ----------------------
  For the Dy (vertical) operators, it is important to basically
  take the diff between two adjacent rows. It is far slower to
  Take a view of a column and compute its diff. In the code below,
  the outer loop walks down the rows, and the inner loop subtracts
  the two adjacent rows. In principle, the inner loop could be replaced
  with something like cblas_daxpy(), but we dont want to overwrite
  our vector.
 */

/*
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
 */
void l1c_Dy(
    l1c_int n, l1c_int m, double alpha, double* restrict X, double* restrict Dy) {
  int row = 0, col = 0;

#pragma omp parallel for private(col)
  for (row = 0; row < n - 1; row++) {
    for (col = row * m; col < (row + 1) * m; col++) {
      Dy[col] = alpha * (X[col + m] - X[col]);
    }
  }

#pragma omp parallel for
  // Ensure the last m items are zero.
  for (int i = (n - 1) * m; i < n * m; i++) {
    Dy[i] = 0;
  }
}

/*
  Given an m by n matrix A, computes lambda*Del_y^T*A. We assume A is stored
  in the 1-D vector A, in row major order. For a 3 x 4 matrix, Del_y^T has
  the matrix representation

 -1     0     0     0     0     0     0     0     0     0     0     0
  0    -1     0     0     0     0     0     0     0     0     0     0
  0     0    -1     0     0     0     0     0     0     0     0     0
  1     0     0    -1     0     0     0     0     0     0     0     0
  0     1     0     0    -1     0     0     0     0     0     0     0
  0     0     1     0     0    -1     0     0     0     0     0     0
  0     0     0     1     0     0    -1     0     0     0     0     0
  0     0     0     0     1     0     0    -1     0     0     0     0
  0     0     0     0     0     1     0     0    -1     0     0     0
  0     0     0     0     0     0     1     0     0     0     0     0
  0     0     0     0     0     0     0     1     0     0     0     0
  0     0     0     0     0     0     0     0     1     0     0     0

 */
void l1c_DyT(
    l1c_int n, l1c_int m, double alpha, double* restrict A, double* restrict dyt) {
  double Ai = 0, Ai_min_m = 0;
  int i = 0, Len = n * m;

#pragma omp parallel for private(Ai, Ai_min_m, i)
  for (i = 0; i < Len; i++) {
    Ai = i < (n - 1) * m ? A[i] : 0;
    Ai_min_m = i < m ? 0 : A[i - m];
    dyt[i] = alpha * (Ai_min_m - Ai);
  }
}

/*
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
 */
void l1c_DyTDy(
    l1c_int n, l1c_int m, double alpha, double* restrict A, double* restrict dytdy) {
  double alp2 = alpha * alpha;
  double Ai_min_m = 0, D_ii_Ai = 0, Ai_p_m = 0;
  int i = 0, Len = n * m;

#pragma omp parallel for private(D_ii_Ai, Ai_min_m, Ai_p_m, i)
  for (i = 0; i < Len; i++) {
    Ai_min_m = i < m ? 0 : A[i - m];
    D_ii_Ai = (i < m || i > m * (n - 1) - 1) ? A[i] : 2 * A[i];
    Ai_p_m = i < Len - m ? A[i + m] : 0;

    dytdy[i] = alp2 * (-Ai_min_m + D_ii_Ai - Ai_p_m);
  }
}
