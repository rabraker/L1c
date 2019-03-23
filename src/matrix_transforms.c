#include "config.h"
#include "l1c_common.h"
#include "l1c_memory.h"

/* Redeclared from l1c_transforms.h*/
typedef struct L1cAxFuns {
  void(*Ax)(double *x, double *y);
  void(*Aty)(double *y, double *x);
  void(*AtAx)(double *x, double *z);

  /* Currently unused.*/
  void(*M)(double *x);
  void(*Mt)(double *y);
  void(*E)(double *x);
  void(*Et)(double *y);
  void(*destroy)(void);

  void *data;
}L1cAxFuns;

/* Forward Declarations */
void Ax(double *x, double *y);
void Aty(double *y, double *x);
void AtAx(double *x_in, double *x_out);
void destroy_matrix_transforms(void);

double *xfm_A=NULL;
l1c_int xfm_N=0;
l1c_int xfm_M=0;
double *xfm_dwork=NULL;


/*
  Do not free A until you are done.
 */
int setup_matrix_transforms(l1c_int n, l1c_int m, double *A, L1cAxFuns *ax_funs){

  int L = n > m ? n : m;

  xfm_dwork = malloc_double(L);
  if (!xfm_dwork){
    return L1C_OUT_OF_MEMORY;
  }
  xfm_A = A;
  xfm_N = n;
  xfm_M = m;

  ax_funs->Ax = Ax;
  ax_funs->Aty = Aty;
  ax_funs->AtAx = AtAx;
  ax_funs->destroy = destroy_matrix_transforms;
  ax_funs->M=NULL;
  ax_funs->Mt=NULL;
  ax_funs->E=NULL;
  ax_funs->Et=NULL;
  ax_funs->data=NULL;

  return 0;
}

void destroy_matrix_transforms(void){
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
  cblas_dgemv(CblasRowMajor, CblasNoTrans, xfm_N, xfm_M, alp, xfm_A, xfm_M, x, inc, beta, y, inc);
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

void AtAx(double *x_in, double *x_out){
  Ax(x_in, xfm_dwork);
  Aty(xfm_dwork, x_out);
}
