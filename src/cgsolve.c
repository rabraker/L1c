/* ##############################################################
This code is a c implementation of the conjugate gradient solver
in l1magic. It is a direct port of the matlab code, which follows
closely to the description on wikipedia

https://en.wikipedia.org/wiki/Conjugate_gradient_method

 Purpose
 -------
 To solve the system of equations
 A x = b
 where A = A^T > 0 via the method of conjugate gradients.
 The advantage to this method is that we do not have to
 store the matrix A, but only need a function which performs
 the linear mapping.


 double *x : result is stored in this array. Should have length n_b.

 double *b : the RHS, should have length n_b.

 size_t n_b: length of b

 double tol : tolerance.

 int verbose : how often to print status updates.

 double *Dwork: A work vector of length 5 * n_b.

 void(*AX_func)(int n, double *x, double *b, void *AX_data) : Pointer to a function
 which evalutes A * x.

 void *AX_data : Data needed by AX_func. For example, if
 you just wanted AX_func to perform a normal matrix multiplication,
 you could do
 AX_func((int n, double *x, double *b, void *AX_data){
 double *A = (double *) AX_data;
 }
*/


#include "cblas.h"
// #include "clapack.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "cgsolve.h"



int cgsolve(double *x, double *b, size_t n_b, double *Dwork,
            void(*AX_func)(int n, double *x, double *b, void *AX_data), void *AX_data,
            CgResults *cg_result, CgParams cg_params){


  int iter, i = 0;

  double delta = 0.0;
  double delta_0= 0.0;
  double delta_old = 0.0;
  double res = 0.0;
  double best_res = 0.0;
  double beta = 0.0;
  double alpha = 0.0;
  double *r, *d, *bestx, *q;

  /* Divide up Dwork for tempory variables */
  r = Dwork;
  d = Dwork + n_b;
  bestx = Dwork + 2 * n_b;
  q = Dwork + 3 * n_b;

  /* Init
  x = zeros(n,1)
  r = b;
  d = r;
  delta = r'*r;
  delta_0 = b'*b;
  numiter = 0;
  bestx = x;
  bestres = sqrt(delta/delta_0);
  */
  for (i=0; i<n_b; i++){
    x[i] = 0.0;
    bestx[i] = 0.0;
  }

  cblas_dcopy((int)n_b, b, 1, r, 1);       /*r=b: copy b (ie, z_i_1) to r */
  cblas_dcopy((int)n_b, r, 1, d, 1);       /*d=r: copy r (ie, z_i_1) to d */
  delta = cblas_ddot(n_b, r, 1, r, 1);     /*delta = r'*r                 */

  delta_0 = cblas_ddot(n_b, b, 1, b, 1);  /*delta_0 = b'*b                */
  best_res = sqrt(delta/delta_0);

  if (cg_params.verbose > 0){
    printf("cg: |Iter| Best resid | Current resid| alpha | beta   |   delta  |\n");
  }
  for (iter=1; iter<=cg_params.max_iter; iter++){

    AX_func(n_b, d, q, AX_data);                 /* q = A * d */

    alpha = delta / cblas_ddot(n_b, d, 1, q, 1); /* alpha delta/(d'*q) */

    cblas_daxpy(n_b, alpha, d, 1, x, 1);         /* x = alpha*d + x    */

    if ( (iter+1 %50 ) == 0){
      AX_func(n_b, x, r, AX_data);               /* r = b - A(x);      */
      cblas_daxpby(n_b, 1.0, b, 1, -1.0, r, 1);  /* r = b - A*x        */
    }else{
      cblas_daxpy(n_b, -alpha, q, 1, r, 1);      /* r = - alpha*q + r; */
    }

    delta_old = delta;
    delta = cblas_ddot(n_b, r, 1, r, 1);         /* delta = r'*r; */

    beta = delta/delta_old;
    cblas_daxpby(n_b, 1.0, r, 1, beta, d, 1);    /* d = r + beta*d; */

    res = sqrt(delta/delta_0);
    if (res < best_res) {
      //bestx = x;
      cblas_dcopy( (int)n_b, x, 1, bestx, 1);
      best_res = res;
    }

    if ( cg_params.verbose >0 && (iter % cg_params.verbose)==0){ // modulo 0 is a floating point exception.
      // printf("cg: Iter = %d, Best residual = %.3e, Current residual = %.3e\n", iter, best_res, res);
      printf("  %d,   %.16e, %.16e, %.16e, %.16e, %.16e  \n", iter, best_res, res, alpha, beta, delta);
    }

    if (delta < (cg_params.tol * cg_params.tol) * delta_0){
      break;
    }

  }


  // x = bestx;
  cblas_dcopy( (int)n_b, bestx, 1, x, 1);
  cg_result->cgres = best_res;
  cg_result->cgiter = iter;


  return 0;

}



/* Simple wrappers around cblas_dgemv to compute b=Ax */

// void dgemv_ColOrder(double *A, int m_A, int n_A, double *x, double *b){
//   enum CBLAS_ORDER MatOrder = CblasColMajor;
//   enum CBLAS_TRANSPOSE NoTrans = CblasNoTrans;
//   cblas_dgemv(MatOrder, NoTrans, m_A, n_A, 1.0, A, m_A, x, 1, 0.0, b, 1);
// }

void dgemv_RowOrder(double *A, int m_A, int n_A, double *x, double *b){
  enum CBLAS_ORDER MatOrder = CblasRowMajor;
  enum CBLAS_TRANSPOSE NoTrans = CblasNoTrans;
  cblas_dgemv(MatOrder, NoTrans, m_A, n_A, 1.0, A, n_A, x, 1, 0.0, b, 1);
}



void Ax_sym(int n, double *x, double *b, void *AX_data){
  double *A = (double *) AX_data;
  enum CBLAS_ORDER Layout = CblasRowMajor;
  enum CBLAS_UPLO uplo = CblasUpper;

  cblas_dspmv (Layout, uplo, n, 1.0, A, x, 1, 0.0, b, 1);

}


void Ax(int n, double *x, double *b, void *AX_data){
  double *A = (double *) AX_data;
  dgemv_RowOrder(A, n, n, x, b);
}


// if (!r){
//   perror("Memory Allocation failure\n");
//   free(r);
//   return 1;
//  }
