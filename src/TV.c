#include "config.h"
#include "l1c_common.h"
#include "l1c_memory.h"
#include <math.h>
#include <stdio.h>

/*
  This file provides a set of functions for computing the discrete
  gradient and laplacian of an N by M image represented as a vector
  in \mathbb{R}^{NM}, which has been concatenated ROW WISE. N.B., that
  this is different from the standard vec(X) operation, which builds a vector
  column wise. However, representing the vector as a row-wise stacking
  conforms to the c convention better.

 */
/* Forward declarations */

void l1c_Dx(l1c_int n, l1c_int m, double *X, double *Dx);

void l1c_DxT(l1c_int n, l1c_int m, double alpha, double *A, double *dxt);

void l1c_DxTDx(l1c_int n, l1c_int m, double alpha, double *A, double *dxtdx);

void l1c_DyT(l1c_int n, l1c_int m, double alpha, double *A, double *dyt);

void l1c_DyTDy(l1c_int n, l1c_int m, double alpha, double *A, double *dytdy);

void l1c_Dy(l1c_int n, l1c_int m, double *X, double *Dy);


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
void l1c_Dx(l1c_int n, l1c_int m, double * restrict X, double * restrict Dx){
  double *X_ = __builtin_assume_aligned(X, DALIGN);
  double *Dx_ = __builtin_assume_aligned(Dx, DALIGN);

  int row=0, col=0;
#pragma omp parallel for private(col)
  for (row=0; row<n; row++){
    for (col=row*m; col<(row+1)*m-1; col++){
      Dx_[col] = X_[col+1] - X_[col];
    }
    // We now have col=(row+1)*m
    Dx_[col] = 0.0;
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
void l1c_DxT(l1c_int n, l1c_int m, double alpha,
             double * restrict A, double * restrict dxt){
  double *A_ = __builtin_assume_aligned(A, DALIGN);
  double *dxt_ = __builtin_assume_aligned(dxt, DALIGN);

  int row, col;
#pragma omp parallel for private(col)
  for (row=0; row<n; row++){
    dxt_[row*m] = -alpha * A_[row*m];
    for (col=row*m+1; col<(row+1)*m-1; col++){
      dxt_[col] = alpha * (A_[col-1] - A_[col]);
    }
    // We now have col=(row+1)*m
    dxt_[col] = alpha * A_[col-1];
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
void l1c_DxTDx(l1c_int n, l1c_int m, double alpha,
               double * restrict A, double * restrict dxtdx){
  double *A_ = __builtin_assume_aligned(A, DALIGN);
  double *dxtdx_ = __builtin_assume_aligned(dxtdx, DALIGN);
  int row, col;

#pragma omp parallel for private(col)
  for (row=0; row<n; row++){
    dxtdx_[row*m] = alpha * (A_[row*m] -A_[row*m+1]);
    for (col=row*m+1; col<(row+1)*m-1; col++){
      dxtdx_[col] = alpha * (-A_[col-1] + 2*A_[col] - A_[col+1]);
    }
    // We now have col=(row+1)*m
    dxtdx_[col] = alpha * (-A_[col-1] + A_[col]);
  }
 }


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
void l1c_Dy(l1c_int n, l1c_int m, double * restrict X, double * restrict Dy){
  double *X_ = __builtin_assume_aligned(X, DALIGN);
  double *Dy_ = __builtin_assume_aligned(Dy, DALIGN);
  int row=0, col=0;

#pragma omp parallel for private(col)
  for (row=0; row<n-1; row++){
    for(col = row * m; col < (row+1) * m; col++){
      Dy_[col] = X_[col+m] - X_[col];
    }
  }

#pragma omp parallel for
  // Ensure the last m items are zero.
  for (int i=(n-1)*m; i<n*m; i++){
    Dy_[i] = 0;
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
void l1c_DyT(l1c_int n, l1c_int m, double alpha,
             double *restrict A, double * restrict dyt){
  double *A_ = __builtin_assume_aligned(A, DALIGN);
  double *dyt_ = __builtin_assume_aligned(dyt, DALIGN);
  double Ai = 0, Ai_min_m = 0;
  int i=0, Len=n*m;

#pragma omp parallel for private(Ai, Ai_min_m, i)
  for (i=0; i<Len; i++){
    Ai = i<(n-1)*m ? A_[i] : 0;
    Ai_min_m = i < m ? 0 : A_[i-m];
    dyt_[i] = alpha * (Ai_min_m - Ai);
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
void l1c_DyTDy(l1c_int n, l1c_int m, double alpha,
               double * restrict A, double * restrict dytdy){
  double *A_ = __builtin_assume_aligned(A, DALIGN);
  double *dytdy_ = __builtin_assume_aligned(dytdy, DALIGN);

  double Ai_min_m=0, D_ii_Ai = 0, Ai_p_m=0;
  int i=0, Len=n*m;

#pragma omp parallel for private(D_ii_Ai, Ai_min_m, Ai_p_m, i)
  for (i=0; i<Len; i++){
    Ai_min_m = i<m ? 0 : A_[i-m];
    D_ii_Ai = (i<m || i > m*(n-1)-1) ? A_[i] : 2 * A_[i];
    Ai_p_m = i < Len-m ? A_[i+m] : 0;

    dytdy_[i] = alpha * (-Ai_min_m + D_ii_Ai - Ai_p_m);
  }

}
