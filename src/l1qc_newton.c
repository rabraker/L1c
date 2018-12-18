#include <stdlib.h>
#include<stdio.h>
#include <math.h>
#include <fftw3.h>
#include "cblas.h"
#include "l1qc_newton.h"
#include "dct.h"
#include "cgsolve.h"

// Need to build these two still
// void A(N, m, z, y);
// void AT(N, m, z, y);

// double fe = 0;
// double *atr, *sigx, *w1p;
#ifndef _MIN_L1QC_
#define _MIN_L1QC_
#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })
#endif

#ifndef _MAX_L1QC_
#define _MAX_L1QC_
#define max(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a > _b ? _a : _b; })
#endif

double sum_abs_vec(int N, double *x){
  int i = 0;
  double sum = 0.0;
  for (i=0; i<N; i++){
    sum = sum + fabs(x[i]);
  }
  return sum;
}


double sum_vec(int N, double *x){
  int i = 0;
  double sum = 0.0;
  for (i=0; i<N; i++){
    sum = sum + x[i];
  }
  return sum;
}

void log_vec(int N, double alpha, double *x, double *logx){
  // Computes alpha*log(x), where x is an array of length N

  int i = 0;
  for (i=0; i<N; i++) {
    logx[i] = log(alpha * x[i] );
  }

}

double logsum(int N, double *x, double alpha) {
  /* Computes sum(log( alpha *x)) */
  int i = 0;
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


/* Computes the H11 part of the hessian. Used to call cgsolve.
 */
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
  double atr_dot_z_fe = 0.0;

  atr_dot_z_fe = cblas_ddot(N, h11p_data.atr, 1, z, 1);
  atr_dot_z_fe = atr_dot_z_fe * h11p_data.one_by_fe_sqrd; //1/fe^2*(atr'*z)

  // y = sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  // y_ = sigx.*z
  for (i=0; i<N; i++){
    y[i] = h11p_data.sigx[i] * z[i];
  }
  // ...+ (1/fe)*At*A*z
  ATA_z = dct_MtEt_EMx_new(z);
  for (i=0; i<N; i++){
    y[i] = y[i] - h11p_data.one_by_fe * ATA_z[i];
  }
  // ...+ 1/fe^2*(atr'*z)*atr;
  for (i=0; i<N; i++){
    y[i] +=  atr_dot_z_fe * h11p_data.atr[i];
  }

}


void get_gradient(int N, double *fu1, double *fu2, double *sigx, double *atr,
                  double fe,  double tau, GradData gd){

  int i = 0;
  double one_by_fu1 = 0.0, one_by_fu2 = 0.0;
  double one_by_fe_Atr = 0.0;
  double ntgu = 0.0, ntgz = 0.0, sig11 = 0.0, sig12 = 0.0;

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

int compute_descent(int N, double *fu1, double *fu2, double *atr, double fe,  double tau,
                    GradData gd, double *Dwork_5N, CgParams cg_params, CgResults *cg_result){
  /* inputs
   --------
   *fu1, *fu2, *r, fe, tau

   cg_params

   outputs
   -------
     dx
     du
     cg_results
     gd.*

  work
  ----
    *Dwork_5N

  */
  int i=0;
  Hess_data h11p_data;
  double *sigx, *Dwork_4N;
  sigx =  Dwork_5N;
  Dwork_4N = Dwork_5N + N;

  get_gradient(N, fu1, fu2, sigx, atr, fe, tau, gd);

  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = 1.0 / (fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;

  cgsolve(gd.dx, gd.w1p, N, Dwork_4N,
          H11pfun, &h11p_data, cg_result, cg_params);

  for (i=0; i<N; i++){
    gd.du[i] = (1.0/gd.sig11[i]) * gd.ntgu[i] - (gd.sig12[i] / gd.sig11[i])  * gd.dx[i];
  }

  return 0;
}

double find_max_step(int N, GradData gd, double *fu1,
                     double *fu2, int M, double *r, double epsilon, int *Iwork_2N){
  double aqe = 0.0, bqe = 0.0, cqe=0.0;
  double smax = 0.0, root = 0.0;
  int *idx_fu1, *idx_fu2;
  int n_idx1 = 0, n_idx2 = 0, idx_i = 0, i=0;
  double min_u1 = 0.0, min_u2=0.0;

  idx_fu1 = Iwork_2N;
  idx_fu2 = Iwork_2N + N;

  aqe = cblas_ddot(M, gd.Adx, 1, gd.Adx, 1);
  bqe = 2.0 * cblas_ddot(M, r, 1, gd.Adx, 1);
  cqe = cblas_ddot(M, r, 1, r, 1) - epsilon * epsilon;

  root = (-bqe + sqrt( bqe*bqe - 4 * aqe * cqe)) / ( 2* aqe);

  /* Find indexes, i.e.,  ifu1 = find((dx-du) > 0);  ifu2 = find((-dx-du) > 0);
   */
  for (i=0; i<N; i++){
    if (gd.dx[i] - gd.du[i] > 0){
      idx_fu1[n_idx1] = i;
      n_idx1++;
    }
  }

  for (i=0; i<N; i++){
    if (-gd.dx[i] - gd.du[i] > 0){
      idx_fu2[n_idx2] = i;
      n_idx2++;
    }
  }


  min_u1 = INFINITY;
  for (i=0; i<n_idx1; i++){
    idx_i = idx_fu1[i];
      min_u1 = min(min_u1, -fu1[idx_i] / (gd.dx[idx_i] - gd.du[idx_i]) );
  }
  min_u2 = INFINITY;
  for (i=0; i<n_idx2; i++){
    idx_i = idx_fu2[i];
      min_u2 = min(min_u2, -fu2[idx_i] / (-gd.dx[idx_i] - gd.du[idx_i]));
  }

  smax = min(min_u1, min_u2);
  smax = min(smax, root);
  return min(1.0, smax) * 0.99;
}


int line_search(int N, int M, double *x, double *u, double *r, double *fu1, double *fu2, GradData gd,
                LSParams ls_params, double *DWORK_5N, double *fe, double *f){
  int iter=0;
  double flin = 0.0, fp = 0.0, fep = 0.0;
  double eps2 = ls_params.epsilon * ls_params.epsilon;
  double a1=0.0, a2 = 0.0, a3 = 0.0;
  double one_by_tau = 1.0/ls_params.tau;
  double step = ls_params.s;
  double rdot = 0.0;
  /* Divy up our work array */
  double *xp = DWORK_5N;
  double *up = DWORK_5N + N;
  double *rp = DWORK_5N + 2*N;
  double *fu1p = DWORK_5N + 3*N;
  double *fu2p = DWORK_5N + 4*N;

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

int newton_init(int N, double *x, double *u,  NewtParams *params,
                double *Dwork,  int M, double *b, int *pix_idx){
  int i = 0;
  double x_max = 0.0, tmp = 0.0;
  double *Ax = Dwork;
  dct_setup(N, M, pix_idx);

  Ax = dct_EMx_new(x);

  cblas_daxpy(M, -1.0, b, 1, Ax, 1); //-b + Ax -->ax
  tmp = cblas_dnrm2(M, Ax, 1)  - params->epsilon;
  if (tmp > 0){
    // Using minimum-2 norm  x0 = At*inv(AAt)*y.') as they
    // will require updates to cgsolve.
    printf("Starting point is infeasible, exiting\n");
    return 1;
  }

  /* Initialize u */
  x_max = -INFINITY;
  for (i=0; i<N; i++){
    x_max = max(x_max, fabs(x[0]));
  }
  for (i=0; i<N; i++){
    u[i] = 0.95 * fabs(x[i]) + 0.10 *x_max;
  }

  /* choose initial value of tau so that the duality gap after the first
     step will be about the original norm */
  //tau = max((2*N+1)/sum(abs(x0)), 1);
  tmp = (double)(2*N+1) / sum_abs_vec(N, x);
  params->tau = max(tmp, 1);

  tmp = log (2 * (double)N +1) - log(params->lbtol) - log(params->tau);
  params->lbiter = ceil (tmp / log(params->mu));

  return 0;
}

int l1qc_newton(int N, double *x, double *u, double *b,
                int M, int *pix_idx, NewtParams params){
  int iter=0, total_newt_iter = 0, tau_iter = 0;
  int status = 0;
  // Line search parameters.
  LSParams ls_params = { .alpha = 0.01,
                         .beta = 0.5,
                         .tau = params.tau,
                         .epsilon = params.epsilon,
                         .s = 1.0};

  double *DWORK_16N = NULL, *DWORK_fftw_2N = NULL;
  int *IWORK_2N = NULL;

  DWORK_16N = calloc(16*N, sizeof(double));
  DWORK_fftw_2N = fftw_alloc_real(2*N);
  IWORK_2N = calloc(2*N, sizeof(int));
  if ( !DWORK_16N | !DWORK_fftw_2N | !IWORK_2N){
    status = 1;
    goto exit;
  }

  double *atr;
  double fe = 0.0, f = 0.0, lambda2 = 0.0, stepsize = 0.0;

  CgParams cg_params = params.cg_params;

  CgResults cg_results;


  double *DWORK_5N = DWORK_16N;
  double *r = DWORK_16N + 5*N;
  double *fu1 = DWORK_16N + 6*N;
  double *fu2 = DWORK_16N + 7*N;

  GradData gd = {.w1p = DWORK_16N   + 8*N,
                 .dx = DWORK_16N    + 9*N,
                 .du = DWORK_16N    + 10*N,
                 .sig11 = DWORK_16N + 11*N,
                 .sig12 = DWORK_16N + 12*N,
                 .ntgu = DWORK_16N  + 13*N,
                 .gradf = DWORK_16N + 14*N,
                 .Adx=NULL}; //gradf needs 2N

  atr = DWORK_fftw_2N;

  /* Compute the initial point. Dwork needs size 2*N */
  f_eval(N, r, x, u, params.tau, params.epsilon, fu1, fu2, &fe, &f, DWORK_5N);

  for (tau_iter=0; tau_iter<params.lbiter; tau_iter++){
    /* ---------------- MAIN ITERATION --------------------- */
    for (iter=0; iter<params.newton_max_iter; iter++){

      /* compute descent direction. returns dx, du, gradf, cgres */
      atr = dct_MtEty(r);
      if(!compute_descent(N, fu1, fu2, atr, fe,  params.tau, gd, DWORK_5N, cg_params, &cg_results)){
        break;
      }

      gd.Adx = dct_EMx_new(gd.dx);

      /* -------------- Line Search --------------------------- */
      ls_params.s = find_max_step(N, gd, fu1, fu2, M, r, params.epsilon, IWORK_2N);

      line_search(N, M, x, u,r, fu1, fu2, gd, ls_params, DWORK_5N, &fe, &f);

      /*Following is for printing only, not calculation. */
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



    /* ----- Update tau or exit ------ */
    total_newt_iter += iter;
    printf("Log barrier iter = %d, l1 = %.3f, functional = %8.3e, tau = %8.3e, total newton-iter =%d\n",
           tau_iter, sum_abs_vec(N, x), sum_vec(N, u), params.tau, total_newt_iter);
    params.tau = params.tau * params.mu;
  }

  goto exit;

  /* ----- Cleanup -------------------- */

 exit:
  free(DWORK_16N);
  fftw_free(DWORK_fftw_2N);
  free(IWORK_2N);
  dct_destroy();

  return status;

}/* MAIN ENDS HERE*/
