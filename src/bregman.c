#include "config.h"
#include <math.h>
#include "l1c_common.h"
#include "bregman.h"
#include "TV.h"
#include "cgsolve.h"
#include <stdio.h>
#include <sys/time.h>
#include <unistd.h>
#include "l1c_math.h"

/* Forward declarations */
static void breg_shrink1(l1c_int N, double *x, double *d, double gamma);

static void breg_mxpy_z(l1c_int N, double * restrict x, double * restrict y, double *z);


/**
   For unit tests, exports a table of function pointers.
 */
BregFuncs breg_get_functions(){
  BregFuncs bfun = {
                    .breg_shrink1 = breg_shrink1,
                    .breg_mxpy_z = breg_mxpy_z};

  return bfun;
}


void breg_shrink1(l1c_int N, double * restrict x, double * restrict d, double gamma){
  double *x_ = __builtin_assume_aligned(x, DALIGN);
  double *d_ = __builtin_assume_aligned(d, DALIGN);
  int i=0;
  double sign_x=0;

#pragma omp parallel for private(sign_x, i)
  for (i=0; i<N; i++){
    sign_x = dsign(x_[i]);
    d_[i] = sign_x * max(fabs(x_[i]) - gamma, 0);
  }
}


/**
   Performs the vector operation
   z = -x + y.
   Assumes x, y, and z are aligned on a DALIGN (default: 64) byte boundary.
*/
void breg_mxpy_z(l1c_int N, double * restrict x, double * restrict y, double *z){
  double *x_ = __builtin_assume_aligned(x, DALIGN);
  double *y_ = __builtin_assume_aligned(y, DALIGN);
  double *z_ = __builtin_assume_aligned(z, DALIGN);

#pragma omp parallel for
  for (int i=0; i<N; i++){
    z_[i] = -x_[i] + y_[i];
  }

}

/**
   Computes the relative 2-norm of error the error between x and y, i.e.,
      nrm = ||x - y||_2/||y||_2
   Assumes x and y are aligned to a DALIGN (default: 64) byte boundary.
 */
double l1c_norm2_err(l1c_int N, double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, DALIGN);
  double *y_ = __builtin_assume_aligned(y, DALIGN);
  int i=0;
  double ynrm = cblas_ddot(N, y, 1, y, 1);
  double nrm = 0.0, diff = 0.0;

  for (i=0; i<N; i++){
    diff = x_[i] - y_[i];
    nrm += diff * diff;
  }
  return sqrt(nrm/ynrm);
}

/*
for n=5, m=4,
D = [2. 3. 3. 2. 3. 4. 4. 3. 3. 4. 4. 3. 3. 4. 4. 3. 2. 3. 3. 2.]
      row
 */
void hess_inv_diag(l1c_int n, l1c_int m, double mu, double lambda, double *D){
  l1c_int k=0, j=0;

  double mu_2lam = 1.0  / (mu + 2*lambda);
  double mu_3lam = 1.0 / (mu + 3*lambda);
  double mu_4lam = 1.0 / (mu + 4*lambda);

  (void)mu;
  (void)lambda;
  // The first m elements.
  D[0] = mu_2lam;
  for(j=1; j<m-1; j++){
    D[k*m + j] = mu_3lam;
  }
  D[m-1] = mu_2lam;

  //The middle (m-2)*n elements
  for(k=1; k<n-1; k++){
    D[k*m] = mu_3lam;
    for(j=1; j<m-1; j++){
      D[k*m + j] = mu_4lam;
    }
    D[k*m+j] = mu_3lam;
  }

  // The last m elements.
  D[(n-1)*m] = mu_2lam;
  for(j=(n-1)*m+1; j<m*n-1; j++){
    D[j] = mu_3lam;
  }
  D[m*n-1] = mu_2lam;
}


void breg_hess_eval(l1c_int N, double *x, double *y, void *hess_data){
  BregTVHessData *BRD =  (BregTVHessData*) hess_data;
  (void)N;
  l1c_int n = BRD->n;
  l1c_int m = BRD->m;

  l1c_DyTDy(n, m, BRD->lambda, x, BRD->dwork1);
  l1c_DxTDx(n, m, BRD->lambda, x, BRD->dwork2);
#pragma omp parallel for
  for(int i=0; i<n*m; i++){
    y[i] = BRD->mu*x[i] + BRD->dwork1[i] + BRD->dwork2[i];
  }
}

void breg_rhs(l1c_int n, l1c_int m, double *f, double *dx, double *bx, double *dy, double *by,
              double *rhs, double mu, double lambda, double *dwork1, double *dwork2){

  l1c_int N = n*m;
  // RHS = mu*f + lam*Delx^T*(dxk-bxk) + lam*Dely^T(dyk - byk)
  breg_mxpy_z(N, by, dy, dwork1);
  l1c_DyT(n, m, lambda, dwork1, dwork2);    //dwork2 = lam*Dely^T(dyk - byk)
  breg_mxpy_z(N, bx, dx, dwork1);
  l1c_DxT(n, m, lambda, dwork1, rhs);       //rhs = lam*Delx^T(xyk - bxk)

  cblas_daxpy(N, mu, f, 1, dwork2, 1);      //dwork2 = mu*f + lam*Dely^T(dyk - byk)
  cblas_daxpy(N, 1.0, dwork2, 1, rhs, 1);   //RHS

}

int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
                       double mu, double tol, l1c_int max_iter){
  long start, end;
  struct timeval timecheck;
  gettimeofday(&timecheck, NULL);
  start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec/1000;


  int iter=0, N=n*m, status=0;
  double lambda = 2*mu;
  double *uk_1=NULL, *d_x=NULL, *d_y = NULL, *b_x=NULL, *b_y=NULL;
  double *dwork1=NULL, *dwork2=NULL, *rhs=NULL, *dwork5N;
  double *Dxu_b=NULL, *Dyu_b=NULL, *D=NULL;
  double dnrm_err = 0;
  double alp;
  uk_1 = malloc_double(N);
  d_x = malloc_double(N);
  d_y = malloc_double(N);
  b_x = malloc_double(N);
  b_y = malloc_double(N);
  Dxu_b = malloc_double(N);
  Dyu_b = malloc_double(N);
  dwork1 = malloc_double(N);
  dwork2 = malloc_double(N);
  rhs = malloc_double(N);
  D = malloc_double(N);
  dwork5N = malloc_double(N*5);
  double *r = dwork5N;
  double *z = dwork5N + N;
  double *q = dwork5N + 2*N;
  if (!uk_1 || !d_x || !d_y || !b_x || !b_y || !Dxu_b || !Dyu_b ||!dwork1 || !dwork2||!rhs ||!D){
    status = 1;
  goto exit;
  }

  hess_inv_diag(n, m, mu, lambda, D);


  BregTVHessData BHD = {.n=n, .m=m, .dwork1=dwork1, .dwork2=dwork2, .mu=mu, .lambda=lambda};
  cblas_dcopy(N, f, 1, uk, 1);
  for (int i=0; i<N; i++){
    uk_1[i] = 0.0;
  }
  for (iter=1; iter<=max_iter; iter++){

    // for (int i=0; i<10; i++){
    breg_rhs(n, m, f, d_x, b_x, d_y, b_y, rhs, mu, lambda, dwork1, dwork2);

    //uk_1 = G^k
    /*Single step of conjugate gradient:*/
    // z0 = M^-1 * r0
    breg_hess_eval(N, uk, r, &BHD);
    cblas_daxpby(N, 1.0, rhs, 1, -1.0, r, 1);       /* r = 1*b + (-1)*Ax           */
    l1c_dxmuly_z(N, r, D, z);
    breg_hess_eval(N, z, q, &BHD);
    alp = cblas_ddot(N, r, 1, z, 1) / cblas_ddot(N, z, 1, q, 1);
    cblas_daxpy(N, alp, z, 1, uk, 1);

    /*Compute Dyu_b = Del_y*u + b */
    l1c_Dx(n, m, uk, Dxu_b);
    cblas_daxpy(N, 1.0, b_x, 1, Dxu_b, 1);
    /*  Dxu_b = Del_x*u + b */
    l1c_Dy(n, m, uk, Dyu_b);
    cblas_daxpy(N, 1.0, b_y, 1, Dyu_b, 1);

    /*Apply shrink operators */
    breg_shrink1(N, Dxu_b, d_x, 1.0/lambda);
    breg_shrink1(N, Dyu_b, d_y, 1.0/lambda);
    // }
    /*Bregman update */
    breg_mxpy_z(N, d_x, Dxu_b, b_x);
    breg_mxpy_z(N, d_y, Dyu_b, b_y);

    dnrm_err = l1c_norm2_err(N, uk, uk_1);
    printf("lambda: %f, dnrm_err: %f\n", lambda, dnrm_err);
    if (dnrm_err < tol){
      break;
    }

    cblas_dcopy(N, uk, 1, uk_1, 1);
  }

  gettimeofday(&timecheck, NULL);
  end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec/1000;
  double time_total = ((double)(end -start)) / 1000.0;
  printf("total c time: %f\n", time_total);

  /*Cleanup before exit. */

 exit:
  free_double(uk_1);
  free_double(d_x);
  free_double(d_y);
  free_double(b_x);
  free_double(b_y);
  free_double(Dxu_b);
  free_double(Dyu_b);

  free_double(dwork1);
  free_double(dwork2);
  free_double(rhs);
  free_double(D);
  free_double(dwork5N);
  return status;
}
