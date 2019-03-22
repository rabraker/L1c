#include "config.h"
#include "l1c_common.h"
#include "l1c_memory.h"


double *xfm_A;
l1c_int xfm_N;
l1c_int xfm_M;
double *xfm_dwork;


/*
  Do not free A until you are done.
 */
int setup_matrix_transforms(l1c_int n, l1c_int m, double *A){

  int L = n > m ? n : m;

  dwork = malloc_double(L);
  if (!dwork){
    return L1C_OUT_OF_MEMORY;
  }
  xfm_A = A;
  xfm_N = n;
  xfm_M = m;

  return 0;
}

void destroy_matrix_transforms(){
  free_double(xfm_dwork);
}

/**
   Computes the matrix-vector product y = A * b, for a full matrix A.
   This is a wrapper for cblas_dgemv.

*/
void Ax(double *x, double *y){
  const double alp = 1.0;
  const double beta = 1.0;
  const l1c_int inc = 1;
  /*For Layout = CblasRowMajor, the value of lda must be at least max(1, m).

    y := alpha*A*x + beta*y

  */
  cblas_dgemv(CblasRowMajor, CblasNoTrans, xfm_N, xfm_M, alp, A, xfm_M, x, inc, beta, y, inc);
}


/**
   Computes the matrix-vector product x = A^T * y, for a full matrix A.

   This is a wrapper for cblas_dgemv.

*/
void Aty(double *y, double *x){
  const double alp = 1.0;
  const double beta = 1.0;
  const l1c_int inc = 1;
  /* dgemv computes
     x = alpha*A^T*y + beta*x, with CblasTrans
   */
  //For Layout = CblasRowMajor, the value of lda must be at least max(1, m).
  cblas_dgemv(CblasRowMajor, CblasTrans, xfm_N, xfm_M, alp, xfm_A, xfm_M, y, inc, beta, x, inc);
}

void AtAx(double x_in, double *x_out){
  const double alp = 1.0;
  const double beta = 1.0;
  const l1c_int inc = 1;

  Ax(x_in, xfm_dwork);
  Aty(xfm_dwork, x_out);
}
