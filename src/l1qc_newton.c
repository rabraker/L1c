#include <stdlib.h>
#include<stdio.h>
#include <math.h>
#include "cblas.h"
#include "l1qc_newton.h"
#include "dct.h"
#include "cgsolve.h"

// Need to build these two still
// void A(N, m, z, y);
// void AT(N, m, z, y);

// double fe = 0;
// double *atr, *sigx, *w1p;



double sum_vec(int N, double *x){
  int i;
  double sum = 0.0;
  for (i=0; i<N; i++){
    sum = sum + x[i];
  }
  return sum;
}

void log_vec(int N, double alpha, double *x, double *logx){
  // Computes alpha*log(x), where x is an array of length N

  int i;
  for (i=0; i<N; i++) {
    logx[i] = log(alpha * x[i] );
  }

}

/* Evalutes the value function */
void f_eval(int N, double *r, double *x, double *u, double tau, double epsilon,
            double *fu1, double *fu2, double *fe, double *f, double *Dwork){
  /*
    inputs
    ------
    N : length of vectors
    r,x,u : input data, of length N
    tau : barrier relaxation paremeter
    epsilon: scalar

    outputs
    -------
    fu1 = x - u
    fu2 = -x - u

    fe  = r^T*r - epsilon ^2
    f = total cost function value

    work space
    ----------
    Dwork : should have length 2*N
   */


  double *log_fu1 = Dwork;
  double *log_fu2 = Dwork+N;

  // fu1 = x - u
  cblas_dcopy(N, x, 1, fu1, 1);
  cblas_daxpy(N, -1.0, u, 1, fu1, 1);

  // fu2 = -x - u
  cblas_dcopy(N, x, 1, fu2, 1);
  cblas_daxpby(N, -1.0, u, 1, -1.0, fu2, 1);

  *fe = 0.5 * (cblas_ddot(N, r, 1, r, 1) - epsilon * epsilon);

  log_vec(N, -1.0, fu1, log_fu1);
  log_vec(N, -1.0, fu2, log_fu2);

  *f = sum_vec(N, u) - (1.0/tau) * (
                               sum_vec(N, log_fu1)
                               + sum_vec(N, log_fu2)
                               + log(-(*fe) )
                                    );

}

// globals sigx, w1p,  atr, fe
// Compute hessian * z

void H11pfun(int N, double *z, double *y,  void *hess_data_in){
  /*
    -- !! Need to have z allocated by fftw_alloc_real(N) !!
    -- Both z and y have dimension N.
    -- Assumes the DCT stuff has been setup already.
   */

  Hess_data hp11_data = *((Hess_data *) hess_data_in);
  // h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  int i = 0;
  double *ATA_z;
  double atr_dot_z_fe;

  atr_dot_z_fe = cblas_ddot(N, hp11_data.atr, 1, z, 1);
  atr_dot_z_fe = atr_dot_z_fe * hp11_data.one_by_fe_sqrd;

  // y = sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  // y_ = sigx.*z
  for (i=0; i<N; i++){
    y[i] = hp11_data.sigx[i] * z[i];
  }
  // ...+ (1/fe)*At*A*z
  ATA_z = dct_MtEt_EMx_new(z);
  for (i=0; i<N; i++){
    y[i] +=  -hp11_data.one_by_fe * ATA_z[i];
  }
  //...+ 1/fe^2*(atr'*z)*atr;
  for (i=0; i<N; i++){
    y[i] +=  atr_dot_z_fe* hp11_data.atr[i];
  }

}


// double flin(double f, double alpha, double s, int N, double *gradf, double *dx, double *du){
//   // For the true linear approximation, alpha = 1.0.
//   // flin = f + alpha*s*(gradf'*[dx; du]);
//   // gradf has length 2*N. dx and du have length N.

//   double delF_z1 = alpha * s * cblas_ddot(N, gradf, 1, dx, 1);
//   double delF_z2 = alpha * s * cblas_ddot(N, gradf+N, 1, du, 1);

//   return f + delF_z1 + delF_z2;
// }


// /*


// globals sigx, w1p,  atr,
void get_gradient(int N, double *fu1, double *fu2, double *wp1, double *sigx, double *atr,
                  double *gradf, double fe,  double tau){

  int i = 0;
  double one_by_fu1, one_by_fu2, one_by_fe_Atr;
  double ntgu, ntgz, sig11, sig12;
  for (i=0; i<N; i++){

    one_by_fu1 = 1.0 / fu1[i];
    one_by_fu2 = 1.0 / fu2[i];
    one_by_fe_Atr = (1.0 / fe ) * atr[i];

    ntgz = one_by_fu1  - one_by_fu2 + one_by_fe_Atr;
    ntgu = -tau - one_by_fu1 - one_by_fu2;

    // gradf = -1/tau*[ntgz; ntgu]
    gradf[i] = (-1.0/tau) * ntgz;
    gradf[i+N] = (-1.0/tau) * ntgu;

    // Hessian parameters
    sig11 = one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;
    sig12 = -one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;

    sigx[i] = sig11 - (sig12 * sig12) / sig11;
    wp1[i] = ntgz - (sig12 * sig12) / sig11;

  }
}




void l1qc_newton(int N, double *x, double *u,
                 void A(int N, int m, double *x, double * y),
                 void At(int N, int m, double *y, double * x),
                 double *b, double epsilon, double tau, double newtontol,
                 int newtonmaxiter, double cgtol, int cgmax_iter, int verbose){

  int iter, line_iter = 0;

  // Line search parameters.
  double alpha = 0.01;
  double beta = 0.5;
  double step_sz = 1.0;
  // x0 = x, u0 = u

  double *r, *fu1, *fu2, *dx, *wp1, *atr, *gradf, *sigx, *Dwork_4n;
  double fe, f = 0.0;
  Hess_data h11p_data;

  Dwork_4n = calloc(4*N, sizeof(double));
  r = calloc(N, sizeof(double));
  fu1 = calloc(N, sizeof(double));
  fu2 = calloc(N, sizeof(double));
  gradf = calloc(2*N, sizeof(double));
  wp1 = calloc(N, sizeof(double));
  dx = calloc(N, sizeof(double));
  sigx = calloc(N, sizeof(double));

  atr = calloc(N, sizeof(double));


  /* Compute the initial point. Dwork needs size 2*N */
  f_eval(N, r, x, u, tau, epsilon, fu1, fu2, &fe, &f, Dwork_4n);

 /* ---------------- MAIN ITERATION --------------------- */
  for (iter=0; iter<newtonmaxiter; iter++){

    atr = dct_MtEty(r);

    get_gradient(N, fu1, fu2, wp1, sigx, atr, gradf, fe, tau);

    h11p_data.one_by_fe = 1.0/fe;
    h11p_data.one_by_fe_sqrd = (1.0/fe) * (1.0/fe);
    h11p_data.atr = atr;
    h11p_data.sigx = sigx;


    cgsolve(dx, wp1, N, cgtol, cgmax_iter, verbose, Dwork_4n,
            H11pfun, &h11p_data);

    /* -------------- Line Search --------------------------- */
    for(line_iter=0; line_iter < 32; line_iter++){

      alpha=alpha*1;
      step_sz = step_sz * beta;
    }

  }




}

