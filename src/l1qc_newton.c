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

double logsum(int N, double *x, double alpha) {
  /* Computes sum(log( alpha *x)) */
  int i=0;
  double total = 0.0;
  for (i=0; i<N; i++){
    total += log(alpha * x[i]);
  }
  return total;
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

  Hess_data h11p_data = *((Hess_data *) hess_data_in);
  // h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  int i = 0;
  double *ATA_z;
  double atr_dot_z_fe;

  atr_dot_z_fe = cblas_ddot(N, h11p_data.atr, 1, z, 1);
  atr_dot_z_fe = atr_dot_z_fe * h11p_data.one_by_fe_sqrd; //1/fe^2*(atr'*z)

  // y = sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  // y_ = sigx.*z
  // for (i=0; i<N; i++){
  //   y[i] = h11p_data.sigx[i] * z[i];
  // }
  // ...+ (1/fe)*At*A*z
  ATA_z = dct_MtEt_EMx_new(z);
  for (i=0; i<N; i++){
    y[i] = ATA_z[i]; //y[i] - h11p_data.one_by_fe * ATA_z[i];
  }
  // //...+ 1/fe^2*(atr'*z)*atr;
  // for (i=0; i<N; i++){
  //   y[i] +=  atr_dot_z_fe * h11p_data.atr[i];
  // }

}



void get_gradient(int N, double *fu1, double *fu2, double *sigx, double *atr,
                  double fe,  double tau, GradData gd){

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
    gd.gradf[i] = (-1.0/tau) * ntgz;
    gd.gradf[i+N] = (-1.0/tau) * ntgu;

    // Hessian parameters
    sig11 = one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;
    sig12 = -one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;

    sigx[i] = sig11 - (sig12 * sig12) / sig11;
    // w1p = ntgz - sig12./sig11.*ntgu;
    gd.w1p[i] = ntgz - (sig12 / sig11) * ntgu;

    gd.sig11[i] = sig11;
    gd.sig12[i] = sig12;
    gd.ntgu[i] = ntgu;
  }
}

int compute_descent(int N, double *fu1, double *fu2, double *r, double fe,  double tau,
                    GradData gd, double *Dwork_5N, CgParams cg_params, CgResults *cg_result){
  int i=0;
  Hess_data h11p_data;
  double *sigx, *Dwork_4N, *atr;
  Dwork_4N = Dwork_5N;
  sigx =  Dwork_5N + 4*N;


  atr = dct_MtEty(r);
  get_gradient(N, fu1, fu2, sigx, atr, fe, tau, gd);

  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = (1.0/fe) * (1.0/fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;

  cgsolve(gd.dx, gd.w1p, N, Dwork_4N,
          H11pfun, &h11p_data, cg_result, cg_params);

  for (i=0; i<N; i++){
  gd.du[i] = (1.0/gd.sig11[i]) * gd.ntgu[i] - (gd.sig12[i] / gd.sig11[i])  * gd.dx[i];

  }
  return 0;
}

double find_max_step(int N, GradData grad_data, double *fu1,
                     double *fu2, double *r, double epsilon){

  return 1.0;
}


int line_search(int N, int M, double *x, double *u, double *r, double *fu1, double *fu2, GradData gd,
                LSParams ls_params, double *DWORK_5N, double *fe, double *f){
  int iter=0;
  double *xp, *up, *rp, *fu1p, *fu2p;
  double flin, fp, fep;
  double eps2 = ls_params.epsilon * ls_params.epsilon;
  double a1, a2, a3;
  double one_by_tau = 1.0/ls_params.tau;
  double step = ls_params.s;
  double rdot;
  /* Divy up our work array */
  xp = DWORK_5N;
  up = DWORK_5N + N;
  rp = DWORK_5N + 2*N;
  fu1p = DWORK_5N + 3*N;
  fu2p = DWORK_5N + 4*N;

  for (iter=0; iter<32; iter++){
    // /* xp = x + s*dx etc*/
    cblas_dcopy(N, x, 1, xp, 1);
    cblas_dcopy(N, u, 1, up, 1);
    cblas_dcopy(M, r, 1, rp, 1);

    cblas_daxpy(N, step, gd.dx, 1, xp, 1);
    cblas_daxpy(N, step, gd.du, 1, up, 1);
    cblas_daxpy(M, step, gd.Adx, 1, rp, 1);

    /* fu1p = xp - up*/
    cblas_dcopy(N, up, 1, fu1p, 1);
    cblas_daxpby(N, 1.0, xp, 1, -1.0, fu1p, 1);

    /* fu2p = -xp - up*/
    cblas_dcopy(N, up, 1, fu2p, 1);
    cblas_daxpby(N, -1.0, xp, 1, -1.0, fu2p, 1);

    rdot = cblas_ddot(M, rp, 1, rp, 1);
    fep = 0.5* (rdot - eps2);

    a1 = logsum(N, fu1p, -1.0);
    a2 = logsum(N, fu2p, -1.0);
    a3 = log(-fep);
    fp = sum_vec(N, up) - one_by_tau *( a1 + a2 +a3);

    /* need gradf'*[dx; du], but dx and du stored separately.
       Use pointer arithmetic with gradf*/

    flin = cblas_ddot(N, gd.gradf, 1, gd.dx, 1);
    flin += cblas_ddot(N, gd.gradf+N, 1, gd.dx, 1);
    flin += *f + ls_params.alpha * step *flin;

    step = ls_params.beta* step;

    printf("iter = %d, fp = %f, flin = %f\n", iter, fp, flin);
    if (fp <= flin){ /* Sufficient decrease test */

      cblas_dcopy(N, xp, 1, x, 1);
      cblas_dcopy(N, up, 1, u, 1);
      cblas_dcopy(M, rp, 1, r, 1);
      cblas_dcopy(N, fu1p, 1, fu1, 1);
      cblas_dcopy(N, fu2p, 1, fu2, 1);
      *f = fp;
      *fe = fep;
      return 0;
    }
  }

  printf("Backtracking line search failed, returning previous iterate.\n");
  printf("Last values: fp = %f, flin =%f\n", fp, flin);

  cblas_dcopy(N, x, 1, xp, 1);
  cblas_dcopy(N, u, 1, up, 1);
  cblas_dcopy(M, r, 1, rp, 1);

  return 1;
}



void l1qc_newton(int N, int M, double *x, double *u,
                 void A(int N, int m, double *x, double * y),
                 void At(int N, int m, double *y, double * x),
                 double *b, NewtParams params){
  // double epsilon, double tau, double newtontol,
  // int newtonmaxiter, double cgtol, int cgmax_iter, int verbose
  int iter;

  // Line search parameters.
  LSParams ls_params;
  ls_params.alpha = 0.01;
  ls_params.beta = 0.5;
  ls_params.tau = params.tau;
  ls_params.epsilon = params.epsilon;
  ls_params.s = 1.0;

  // double cgres = 0.0;
  // x0 = x, u0 = u

  double *r, *fu1, *fu2, *DWORK_5N;
  double fe, f, lambda2, stepsize = 0.0;

  GradData gd;
  CgParams cg_params;
  CgResults cg_results;

  DWORK_5N = calloc(5*N, sizeof(double));
  r = calloc(N, sizeof(double));
  fu1 = calloc(N, sizeof(double));
  fu2 = calloc(N, sizeof(double));

  gd.w1p = calloc(N, sizeof(double));
  gd.dx = calloc(N, sizeof(double));
  gd.du = calloc(N, sizeof(double));
  gd.gradf = calloc(2*N, sizeof(double));
  gd.Adx = calloc(2*N, sizeof(double));
  gd.sig11 = calloc(N, sizeof(double));
  gd.sig12 = calloc(N, sizeof(double));
  gd.w1p = calloc(N, sizeof(double));
  gd.ntgu = calloc(N, sizeof(double));




  /* Compute the initial point. Dwork needs size 2*N */
  f_eval(N, r, x, u, params.tau, params.epsilon, fu1, fu2, &fe, &f, DWORK_5N);

 /* ---------------- MAIN ITERATION --------------------- */
  for (iter=0; iter<params.newton_max_iter; iter++){

    /* compute descent direction. returns dx, du, gradf, cgres */
    if(!compute_descent(N, fu1, fu2, r, fe,  params.tau, gd, DWORK_5N, cg_params, &cg_results)){
      return;
    }


    ls_params.s = find_max_step(N, gd, fu1, fu2, r, params.epsilon);

    /* -------------- Line Search --------------------------- */
    line_search(N, M, x, u,r, fu1, fu2, gd, ls_params, DWORK_5N, &fe, &f);

    lambda2 = cblas_ddot(N, gd.gradf, 1, gd.dx, 1);
    lambda2 += cblas_ddot(N, gd.gradf+N, 1, gd.du, 1);
    /* want norm [dx; du] */
    stepsize = cblas_ddot(N, gd.dx, 1, gd.dx, 1);
    stepsize += cblas_ddot(N, gd.du, 1, gd.du, 1);
    stepsize = sqrt(stepsize) * 1.0; //* (ls_params.s);

    if (params.verbose >0){
      printf("Newton iter = %d, Functional = %8.3f, Stepsize = %8.3f",
             iter, -lambda2/2.0, stepsize);
      printf("                  CG RES = %8.3f, CG iter = %d\n",
             cg_results.cgres, cg_results.cgiter);

    }

  }/* main iter*/

}/* MAIN ENDS HERE*/
