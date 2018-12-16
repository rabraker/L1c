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

 double *Dwork: A work vector of length 4 * n_b.

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
  //r=b
  cblas_dcopy((int)n_b, b, 1, r, 1);    /* copy b (ie, z_i_1) to r */
  //d=r
  cblas_dcopy((int)n_b, r, 1, d, 1);    /* copy r (ie, z_i_1) to d */
  // cblas_dcopy((int)n_b, r, 1, q, 1);    /* copy r (ie, z_i_1) to d */
  //delta = r'*r
  delta = cblas_ddot(n_b, r, 1, r, 1);

  //delta_0 = b'*b
  delta_0 = cblas_ddot(n_b, b, 1, b, 1);
  best_res = sqrt(delta/delta_0);


  for (iter=0; iter<cg_params.max_iter; iter++){
    // void catlas_daxpby(const int __N, const double __alpha, const double *__X, const int __incX,
    // const double __beta, double *__Y, const int __incY);

    AX_func(n_b, d, q, AX_data);

    // alpha delta/(d'*q);
    alpha = delta / cblas_ddot(n_b, d, 1, q, 1);

    // x = alpha*d + x;
    cblas_daxpy(n_b, alpha, d, 1, x, 1);

    if ( (iter+1 %50 ) == 0){
      // r = b - A(x);
      AX_func(n_b, x, r, AX_data);
      // r = b - A*x
      cblas_daxpby(n_b, 1.0, b, 1, -1.0, r, 1);
    }else{
      // r = - alpha*q + r;
      cblas_daxpy(n_b, -alpha, q, 1, r, 1);
    }

    delta_old = delta;
    //delta = r'*r;
    delta = cblas_ddot(n_b, r, 1, r, 1);

    beta = delta/delta_old;
    // d = r + beta*d;
    cblas_daxpby(n_b, 1.0, r, 1, beta, d, 1);

    res = sqrt(delta/delta_0);
    if (res < best_res) {
      //bestx = x;
      cblas_dcopy( (int)n_b, x, 1, bestx, 1);
      best_res = res;
    }

    if ( cg_params.verbose >0 && (iter % cg_params.verbose)==0){ // modulo 0 is a floating point exception.
      printf("cg: Iter = %d, Best residual = %.3f, Current residual = %.3f\n",
             iter, best_res, res);
    }

    if (delta < pow(cg_params.tol, 2) * delta_0){
      break;
    }

  }


  // x = bestx;
  cblas_dcopy( (int)n_b, bestx, 1, x, 1);
  cg_result->cgres = best_res;
  cg_result->cgiter = iter;

  printf("CG-iters = %d\n", iter);

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



void Ax(int n, double *x, double *b, void *AX_data){
  double *A = (double *) AX_data;
  dgemv_RowOrder(A, n, n, x, b);
}


// if (!r){
//   perror("Memory Allocation failure\n");
//   free(r);
//   return 1;
//  }
