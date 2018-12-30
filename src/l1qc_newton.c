#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#include "l1qc_newton.h"
#include "dct.h"
#include "cgsolve.h"
#include "l1qc_common.h"



void axpy_z(size_t N, double alpha, double * restrict x, double * restrict y, double * restrict z){
  /* Computes z = a * x + y. Similary to cblas_axpy, but for when you don't want to overwrite y.
   This way, we avoid a call to cblas_dcopy().

   May be worth explicitly vectorizing (e.g., ispc??) this function.
  */
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);
  double *z_ = __builtin_assume_aligned(z, 64);

  size_t i;
  for (i = 0; i<N; i++){
    z_[i] = alpha * x_[i] + y_[i];
  }
}


double sum_abs_vec(int N, double *x){
  return cblas_dasum(N, x, 1);
}


double sum_vec(int N, double *x){
  int i = 0;
  double sum = 0.0;
  for (i=0; i<N; i++){
    sum = sum + x[i];
  }
  return sum;
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
void f_eval(int N, double *x, double *u, int M, double *r, double tau, double epsilon,
            double *fu1, double *fu2, double *fe, double *f){
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

  double a1=0, a2=0, a3=0;
  // double *log_fu1 = Dwork_2N;
  // double *log_fu2 = Dwork_2N + N;

  // fu1 = x - u
  cblas_dcopy(N, x, 1, fu1, 1);
  cblas_daxpy(N, -1.0, u, 1, fu1, 1);

  // fu2 = -x - u
  cblas_dcopy(N, x, 1, fu2, 1);
  cblas_daxpby(N, -1.0, u, 1, -1.0, fu2, 1);

  *fe = 0.5 * (cblas_ddot(M, r, 1, r, 1) - epsilon * epsilon);

  a1 = logsum(N, fu1, -1.0);
  a2 = logsum(N, fu2, -1.0);
  a3 = log(- (*fe));
  *f = sum_vec(N, u) - (1.0/tau) * ( a1 + a2 +a3);

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
                     double *fu2, int M, double *r, double epsilon){
  double aqe = 0.0, bqe = 0.0, cqe=0.0;
  double smax = 0.0, root = 0.0;
  size_t i=0;
  double min_u1 = 0.0, min_u2=0.0;

  aqe = cblas_ddot(M, gd.Adx, 1, gd.Adx, 1);
  bqe = 2.0 * cblas_ddot(M, r, 1, gd.Adx, 1);
  cqe = cblas_ddot(M, r, 1, r, 1) - epsilon * epsilon;

  root = (-bqe + sqrt( bqe*bqe - 4 * aqe * cqe)) / ( 2* aqe);

  /*Towards the end, the root becomes complex. When that happens, cant take
   the real part, because that is just -bqe which is less than zero.*/
  if ( isnan(root) ){
    PRINT("Warning: maximum step which satisfies cone constraints has become complex.\n");
    PRINT("Trying smax=1.0 for the cone-constraing portion.\n");
    root = 1.0;
  }

  /* Implements the matlab code:
     ifu1 = find((dx-du) > 0);
     ifu2 = find((-dx-du) > 0);
     min_u1 =  -fu1(ifu1)./(dx(ifu1)-du(ifu1));
     min_u2 = -fu2(ifu2)./(-dx(ifu2)-du(ifu2)); ...

  */
  min_u1 = INFINITY;
  for (i=0; i<N; i++){
    if ( !(gd.dx[i] - gd.du[i] > 0) ){
      continue;
    }else{
    min_u1 = min(min_u1, -fu1[i] / (gd.dx[i] - gd.du[i]) );
    }
  }

  min_u2 = INFINITY;
  for (i=0; i<N; i++){
    if ( !(-gd.dx[i] - gd.du[i] > 0) ){
      continue;
    }else{
      min_u2 = min(min_u2, -fu2[i] / (-gd.dx[i] - gd.du[i]));
    }
  }



  smax = min(min_u1, min_u2);
  smax = min(smax, root);
  return min(1.0, smax) * 0.99;
}


LSStat line_search(int N, int M, double *x, double *u, double *r, double *b, double *fu1, double *fu2, GradData gd,
                LSParams ls_params, double *DWORK_5N, double *fe, double *f){
  LSStat ls_stat = {.flx=0, .flu = 0, .flin=0, .step=0, .status=0};
  int iter=0;
  double flx=0, flu=0, flin = 0.0, fp = 0.0, fep = 0.0;

  double *fftw_tmp;
  double step = ls_params.s;
  //double rdot = 0.0;
  /* Divy up our work array */
  double *xp = DWORK_5N;
  double *up = DWORK_5N + N;
  double *rp = DWORK_5N + 2*N;
  double *fu1p = DWORK_5N + 3*N;
  double *fu2p = DWORK_5N + 4*N;

  for (iter=1; iter<=32; iter++){
    // /* xp = x + s*dx etc*/
    axpy_z(N, step, gd.dx, x, xp);
    axpy_z(N, step, gd.du, u, up);

    /*Re-evaluate, rather than relying on Adx */
    fftw_tmp = dct_EMx_new(xp);
    axpy_z(M, -1.0, b, fftw_tmp, rp); //rp = A*xp - b

    f_eval(N, xp, up, M, rp, ls_params.tau, ls_params.epsilon, fu1p, fu2p, &fep, &fp);

    /* need gradf'*[dx; du], but dx and du stored separately.
       Use pointer arithmetic with gradf*/

    flx = cblas_ddot(N, gd.gradf, 1, gd.dx, 1);
    flu = cblas_ddot(N, gd.gradf+N, 1, gd.du, 1);
    flin = *f + ls_params.alpha * step * (flx + flu);

    //PRINT("iter = %d, fp = %f, flin = %f\n", iter, fp, flin);
    if (fp <= flin){ /* Sufficient decrease test */

      cblas_dcopy(N, xp, 1, x, 1);
      cblas_dcopy(N, up, 1, u, 1);
      cblas_dcopy(M, rp, 1, r, 1);
      cblas_dcopy(N, fu1p, 1, fu1, 1);
      cblas_dcopy(N, fu2p, 1, fu2, 1);
      *f = fp;
      *fe = fep;
      ls_stat.flx = flx;
      ls_stat.flu = flu;
      ls_stat.flin = flin;
      ls_stat.step = step;
      ls_stat.iter = iter;
      ls_stat.status = 0;
      return ls_stat;
    }

    step = ls_params.beta * step;
  }

  PRINT("Backtracking line search failed, returning previous iterate.\n");
  PRINT("Last values: fp = %f, flin =%f\n", fp, flin);

  cblas_dcopy(N, x, 1, xp, 1);
  cblas_dcopy(N, u, 1, up, 1);
  cblas_dcopy(M, r, 1, rp, 1);

  ls_stat.flx = flx;
  ls_stat.flu = flu;
  ls_stat.flin = flin;
  ls_stat.step = step;
  ls_stat.iter = min(iter, 32);
  ls_stat.status = 1;
  return ls_stat;
}

int newton_init(int N, double *x, double *u,  NewtParams *params,
                int M, int *pix_idx){
  size_t i = 0;
  double x_max = 0.0, tmp = 0.0;

  dct_setup(N, M, pix_idx);

  /* Initialize u */
  x_max = -INFINITY;
  for (i=0; i<N; i++){
    x_max = max(x_max, fabs(x[i]));
  }
  for (i=0; i<N; i++){
    u[i] = 0.95 * fabs(x[i]) + 0.10 *x_max;
  }

  /* choose initial value of tau so that the duality gap after the first
     step will be about the original norm */
  //tau = max((2*N+1)/sum(abs(x0)), 1);
  tmp = sum_abs_vec(N, x);
  tmp = (double)(2*N+1) / tmp;
  params->tau = max(tmp, 1);

  tmp = log (2 * (double)N + 1) - log(params->lbtol) - log(params->tau);
  if (params->lbiter == 0) // If already set, use that instead.
    params->lbiter =(int) ceil (tmp / log(params->mu));

  return 0;
}

int check_feasible_start(int M, double *r, double epsilon){
  double tmp = cblas_dnrm2(M, r, 1)  - epsilon;
  if (tmp > 0){
    // Using minimum-2 norm  x0 = At*inv(AAt)*y.') as they
    // will require updates to cgsolve.
    PRINT("epsilon = %f,  ||r|| = %.20f\n", epsilon, cblas_dnrm2(M, r, 1));
    return 1;
  }else
    return 0;
}


int save_x(int N, double *x, char *fname){
  FILE *fid = fopen(fname, "w");
  if (fid == NULL){
    PRINT("Error opening %s\n", fname);
    return 1;
  }

  for (int i=0; i<N-1; i++){
    fprintf(fid, "%.15f, ", x[i]);
  }
  fprintf(fid, "%.15f ", x[N-1]);
  return 0;
}

int l1qc_newton(int N, double *x, double *u, double *b,
                int M, int *pix_idx, NewtParams params){
  LSStat ls_stat;// = {.flx=0, .flu = 0, .flin=0, .step=0, .status=0};
  CgParams cg_params = params.cg_params;
  CgResults cg_results;

  int iter=0, total_newt_iter = 0, tau_iter = 0;
  int status = 0;

  double fe = 0.0, f = 0.0, lambda2 = 0.0, stepsize = 0.0;

  double *atr=NULL, *DWORK_16N = NULL, *DWORK_fftw_2N = NULL;
  double *fftw_tmp;
  DWORK_16N = malloc_double(16*N);
  DWORK_fftw_2N = fftw_alloc_real(2*N);

  if ( !DWORK_16N | !DWORK_fftw_2N){
    status = 1;
    goto exit;
  }

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

  gd.Adx = malloc_double(M);

  atr = DWORK_fftw_2N+N;

  newton_init(N, x, u, &params, M, pix_idx);

  // Line search parameters.
  LSParams ls_params = { .alpha = 0.01,
                         .beta = 0.5,
                         .tau = params.tau,
                         .epsilon = params.epsilon,
                         .s = 1.0};

  if (params.verbose >0) {
    PRINT("Total Log-Barrier iterations:  %d \n", params.lbiter);
  }

  /* ---------------- MAIN **TAU** ITERATION --------------------- */
  for (tau_iter=1; tau_iter<=params.lbiter; tau_iter++){
    if (params.verbose > 1){
      PRINT("Newton-iter | Functional | Newton decrement |  Stepsize  |  cg-res | cg-iter | backiter |    s    | \n");
      PRINT("----------------------------------------------------------------------------------------------------\n");
    }
    /* Compute Ax - b = r */
    fftw_tmp = dct_EMx_new(x);
    cblas_dcopy(M, fftw_tmp, 1, r, 1);
    cblas_daxpy(M, -1.0, b, 1, r, 1); //-b + Ax -->ax

    if ( (tau_iter==1) & check_feasible_start(M, r, params.epsilon) ){
        PRINT("Starting point is infeasible, exiting\n");
        status = 1;
        goto exit;
      }

    f_eval(N, x, u, M, r, params.tau, params.epsilon, fu1, fu2, &fe, &f);

    /* ---------------- MAIN Newton ITERATION --------------------- */
    for (iter=1; iter<+params.newton_max_iter; iter++){

     /* Note: we do not need to re-evaluate the functional because we now do
        that inside the line search. Thus, f, fe, fu1, fu2 and r provided by
        linesearch should be exact given the current stepsize.
      */

      /* compute descent direction. returns dx, du, gradf, cgres */
      fftw_tmp = dct_MtEty(r);
      cblas_dcopy(N, fftw_tmp, 1, atr, 1);

      if(compute_descent(N, fu1, fu2, atr, fe,  params.tau, gd, DWORK_5N, cg_params, &cg_results)){
        break;
      }

      fftw_tmp= dct_EMx_new(gd.dx);
      cblas_dcopy(M, fftw_tmp, 1, gd.Adx, 1);
      /* -------------- Line Search --------------------------- */
      ls_params.s = find_max_step(N, gd, fu1, fu2, M, r, params.epsilon);

      ls_stat = line_search(N, M, x, u,r,b, fu1, fu2, gd, ls_params, DWORK_5N, &fe, &f);
      if (ls_stat.status > 0)
        break;


      /* Check for exit */
      lambda2 = -(cblas_ddot(N, gd.gradf, 1, gd.dx, 1) + cblas_ddot(N, gd.gradf+N, 1, gd.du, 1) );

      if (params.verbose >1){
        /*Following is for printing only, not calculation. */

        /* want norm [dx; du] */
        stepsize = cblas_ddot(N, gd.dx, 1, gd.dx, 1) + cblas_ddot(N, gd.du, 1, gd.du, 1);
        stepsize = stepsize * ls_params.s;

        /*            NI         fcnl         dec         sz     cgr       cgI        BI       s  */
        PRINT("     %3d       %8.3e       %08.3e    % 8.3e   %08.3e     %3d       %2d      %.3g \n",
                     iter, f, lambda2/2, stepsize, cg_results.cgres, cg_results.cgiter, ls_stat.iter,
               ls_stat.step);
#ifdef __MATLAB__
        mexEvalString("drawnow('update');");
#endif
      }

      /*Check for early exit */
      if (lambda2/2 < params.newton_tol){
        break;
      }


    }/* main iter*/



    /* ----- Update tau or exit ------ */
    total_newt_iter += iter;
    if (params.verbose > 0){
      PRINT("\n********************************************************************************************\n");
      PRINT("Log barrier iter = %d, l1 = %.3f, functional = %8.3e, tau = %8.3e, total newton-iter =%d\n",
            tau_iter, sum_abs_vec(N, x), sum_vec(N, u), params.tau, total_newt_iter);
      PRINT("********************************************************************************************\n\n");
    }
#ifdef __MATLAB__
    mexEvalString("drawnow('update');");
#endif
    params.tau = params.tau * params.mu;
    ls_params.tau = params.tau;
  } /* LB-iter */

  /* ----- Cleanup -------------------- */
  goto exit;

 exit:
  free_double(DWORK_16N);
  free_double(gd.Adx);
  dct_destroy();

  return status;

}/* MAIN ENDS HERE*/
