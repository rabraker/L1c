/**  \file
This contains the conjugate gradient solver, cgsolve. The two small routines Ax and Ax_sym illustrate how the user function AX_func can parse the input void *AX_data.

*/

#ifdef _USEMKL_
#include "mkl.h"
#else
#include "cblas.h"
#endif
// #include "clapack.h"
#include <stdlib.h>
#include <math.h>
// #include <stdio.h>
#include "cgsolve.h"
#include "l1c_common.h"

/**

   Purpose
   -------
   To solve the system of equations
   \f[
   A x = b
   \f]
   where \f$A = A^T > 0\f$ via the method of conjugate gradients.
   The advantage to this method is that we do not have to
   store the matrix A, but only need a function which performs
   the linear mapping.


   @param[out] x The result is stored in this array. Should have length n_b.
   @param[in] b  The RHS, should have length n_b.
   @param[in] n_b Length of the vector b.
   @param[in] Dwork Pointer to a work array of length 5 * n_b.
   @param[in] void(*AX_func) Pointer to a function which evalutes A * x.
   @param[in] AX_data  Data needed by AX_func. For example, if
   you just wanted AX_func to perform a normal matrix multiplication,
   you could do
   AX_func((int n, double *x, double *b, void *AX_data){
   double *A = (double *) AX_data;
   @param[out] cg_result Pointer to a struct containing relevent results from the computation.
   @param[in] cg_params Struct containing parameters (tolerance and verbosity) for the computation.


   Algorithm
   ---------
   This code is a c implementation of the conjugate gradient solver
   in l1magic. It is a direct port of the matlab code, which follows
   closely to the description on wikipedia

   https://en.wikipedia.org/wiki/Conjugate_gradient_method

 }
*/

int cgsolve(double *x, double *b, l1c_int n_b, double *Dwork,
            void(*AX_func)(l1c_int n, double *x, double *b, void *AX_data), void *AX_data,
            CgResults *cg_result, CgParams cg_params){


  int iter;
  l1c_int i = 0;

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
  cg_result->cgiter = min(iter, cg_params.max_iter); //Loops increment before exiting.


  return 0;

}



/**
   Computes the matrix-vector product y = A * b, for a symmetric matrix A.
   This is a wrapper for cblas_dspmv.
*/

void Ax_sym(l1c_int n, double *x, double *b, void *AX_data){

  double *A = (double *) AX_data;

  cblas_dspmv (CblasRowMajor, CblasUpper, n, 1.0, A, x, 1, 0.0, b, 1);

}

/**
Computes the matrix-vector product y = A * b, for a full matrix A.
This is a wrapper for cblas_dgemv.
 */
void Ax(l1c_int n, double *x, double *b, void *AX_data){
  double *A = (double *) AX_data;

  cblas_dgemv(CblasRowMajor, CblasNoTrans, n, n, 1.0, A, n, x, 1, 0.0, b, 1);
}
