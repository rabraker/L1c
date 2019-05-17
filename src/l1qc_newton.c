#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cblas.h"
#include "l1c.h"
#include "vcl_math.h"
#include "l1c_math.h"
#include "l1qc_newton.h"


/* ---------------- Forward Declarations ---------------------- */

static void lb_report(int lb_iter, int m, double *u, double l1,
                      l1c_L1qcOpts params, l1c_LBResult lb_res);

static void newton_report(int iter, int m, l1c_l1qcProb *Prb, double lambda2,
                          l1c_CgResults cg_results, LSStat ls_status);



/* Return pointers to all the static functions. To by called by testing code only.*/


/* Evalutes the value functional */
void _l1c_l1qc_f_eval( l1c_l1qcProb *Prb, double *x, double *u){
  /*
    inputs
    ------
    m : length of vectors
    r,x,u : input data, of length m
    tau : barrier relaxation paremeter
    epsilon: scalar

    outputs
    -------
    fu1 = x - u
    fu2 = -x - u

    fe  = <r, r> - epsilon ^2
    f = total cost function value

    work space
    ----------
    Dwork : should have length 2*m
   */

  double a1=0, a2=0, a3=0;

  double epsilon = Prb->epsilon;
  l1c_int m = Prb->m;
  l1c_int n = Prb->n;

  /* Compute Ax - b = r */
  Prb->ax_funs.Ax(x, Prb->r);
  cblas_daxpy(n, -1.0, Prb->b, 1, Prb->r, 1);

  // fu1 = x - u
  l1c_daxpy_z(m, -1.0, u, x, Prb->fu1);
  // fu2 = -x - u
  axpby_z(m, -1.0, x, -1.0, u, Prb->fu2);

  Prb->fe_val = 0.5 * (cblas_ddot(n, Prb->r, 1, Prb->r, 1) - epsilon * epsilon);

  a1 = vcl_logsum(m, -1.0, Prb->fu1);
  a2 = vcl_logsum(m, -1.0, Prb->fu2);
  a3 = log(- (Prb->fe_val));
  Prb->f_val = vcl_sum(m, u) - (1.0/Prb->tau) * ( a1 + a2 +a3);

}


/* Computes the H11 part of the hessian. Used to call cgsolve.
 */
void _l1c_l1qc_H11pfun(l1c_int m, double *z, double *y,  void *hess_data_in){
  /*

   */

  // Hess_data h11p_data = *((Hess_data *) hess_data_in);
  // // // h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;

  // double *ATA_z = h11p_data.Dwork_1m;
  // double atr_dot_z_fe = 0.0;

  // atr_dot_z_fe = cblas_ddot(m, h11p_data.atr, 1, z, 1);
  // atr_dot_z_fe = atr_dot_z_fe * h11p_data.one_by_fe_sqrd; //1/fe^2*(atr'*z)


  // h11p_data.AtAx(z, ATA_z);                               // AtA_z
  // l1c_dxmuly_z(m, h11p_data.sigx, z, y);                  // y = sigx.*z
  // cblas_daxpy(m, -h11p_data.one_by_fe, ATA_z, 1, y, 1);   // y = sigx.*z (1/fe)*ATA_z
  // cblas_daxpy(m, atr_dot_z_fe, h11p_data.atr, 1, y, 1);

  Hess_data h11p_data = *((Hess_data *) hess_data_in);
  // h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;

  double atr_dot_z_fe = 0.0;

  atr_dot_z_fe = cblas_ddot(m, h11p_data.atr, 1, z, 1);
  atr_dot_z_fe = atr_dot_z_fe * h11p_data.one_by_fe_sqrd; //1/fe^2*(atr'*z)


  // h11pfun = @(z) sigx.*z - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
  h11p_data.AtAx(z, y);                               // y = A^T*A*z
  l1c_daxpby(m, atr_dot_z_fe, h11p_data.atr,          // y = - (1/fe)*At(A(z)) + 1/fe^2*(atr'*z)*atr;
               -h11p_data.one_by_fe, y);
  vcl_dxMy_pz(m, h11p_data.sigx, z, y);                // y = sigx.*z + A^T*A*z


}


/*

For f1 = x-u
    f2 = -x-u
    f3 = ||Ax-b||^2 - eps^2

and cost f(x,u) = (0,1)^T(x,u) + (1/tau) sum( -log(-f1))
                  +(1/tau) sum(-log(-f2))
                  +(1/tau) sum(-log(f3))

the negative gradient multiplied by tau is partitioned as (delf_x, delf_u) where

    neg_tau_delx = [ 1/f1  - 1/f2 + (1/f3)*A^Tr     (where r = Ax-b)
    neg_tau_delu = [ 1 - 1/f1 - 1/f2]

The hessian is given by Hx + Hu where

         |-(1/f3)A^TA + (1/f3^2)A^Trr^TA  |   0   |
    Hx = |--------------------------------|-------|
         |              0                 |   0   |

    Hu = | (1/f1)^2 + (1/f2)^2  |

    |                | diag(-1/f1^2) |
    |diag(1/f1^2)    | + diag(1/f2^2)|
    |  + diag(1/f2^2)|               |   |Hu_11 | Hu_12 |
    |--------------------------------| = |--------------|
    | diag(-1/f1^2)  |diag(1/f1^2)   |   |Hu_12 | Hu_11 |
    |  +diag(1/f2^2) | + diag(1/f2^2)|
    |                |               |

 */
void _l1c_l1qc_hess_grad(l1c_l1qcProb *Prb, double *sigx, double *atr){

  l1c_int i = 0;
  double one_by_fu1 = 0.0, one_by_fu2 = 0.0;
  double one_by_fe_Atr = 0.0;
  double ntgu = 0.0, ntgz = 0.0, sig11 = 0.0, sig12 = 0.0;
  double mone_by_tau = -1.0/Prb->tau;
  double tau = Prb->tau;
  double fe = Prb->fe_val;

  #pragma omp parallel for private(one_by_fu1, one_by_fu2, one_by_fe_Atr, ntgz, ntgu, sig11, sig12)
  for (i=0; i < Prb->m; i++){

    one_by_fu1 = 1.0 / Prb->fu1[i];
    one_by_fu2 = 1.0 / Prb->fu2[i];
    one_by_fe_Atr = (1.0 / fe ) * atr[i];

    ntgz = one_by_fu1  - one_by_fu2 + one_by_fe_Atr;
    ntgu = -tau - one_by_fu1 - one_by_fu2;

    // gradf = -1/tau*[ntgz; ntgu]
    Prb->gradf[i] = mone_by_tau * ntgz;
    Prb->gradf[i + Prb->m] = mone_by_tau * ntgu;

    // Hessian parameters
    sig11 = one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;
    sig12 = -one_by_fu1 * one_by_fu1 + one_by_fu2 * one_by_fu2;

    sigx[i] = sig11 - (sig12 * sig12) / sig11;
    // w1p = ntgz - sig12./sig11.*ntgu;
    Prb->w1p[i] = ntgz - (sig12 / sig11) * ntgu;

    Prb->sig11[i] = sig11;
    Prb->sig12[i] = sig12;
    Prb->ntgu[i] = ntgu;
  }
}


/*
  Computes the descent direction (dx,du) as the solution to
  H (dx,du) = -(\nabla_x f(x,u), \nabla_u f(x,u)), where H is the Hessian. We decompose
  the hessian.

  The Hessian is split into 4 blocks: only the H11 part is dense.
  The rest (H12, H21, and H22) are diagonal.
 */
int _l1c_l1qc_descent_dir(l1c_l1qcProb *Prb, l1c_CgParams cg_params, l1c_CgResults *cg_result){
  /* inputs
   --------
   *fu1, *fu2, *atr, fe, tau

   cg_params

   outputs
   -------
     dx
     du
     cg_results
     gd.*

  work
  ----
    *Dwork_5m

  */
  l1c_int i=0;
  Hess_data h11p_data;
  double *sigx, *atr;
  double **Dwork_4m;
  double fe = Prb->fe_val;

  sigx =  Prb->DWORK7[0];
  h11p_data.Dwork_1m = Prb->DWORK7[1];
  atr = Prb->DWORK7[2];
  Dwork_4m = Prb->DWORK7 + 3;

  Prb->ax_funs.Aty(Prb->r, atr); //atr = A'*r.

  _l1c_l1qc_hess_grad(Prb, sigx, atr);

  h11p_data.one_by_fe = 1.0 / fe;
  h11p_data.one_by_fe_sqrd = 1.0 / (fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;
  h11p_data.AtAx = Prb->ax_funs.AtAx;

  l1c_cgsolve(Prb->m, Prb->dx, Prb->w1p, Dwork_4m,
              _l1c_l1qc_H11pfun, &h11p_data, cg_result, cg_params);

  if (cg_result->cgres > 0.5){
    return L1C_CGSOLVE_FAILURE;
  }

#pragma omp parallel for
  for (i=0; i < Prb->m; i++){
    // gd.du[i] = (gd.ntgu[i] - gd.sig12[i] * gd.dx[i] ) / gd.sig11[i];
    Prb->du[i] = (1.0/Prb->sig11[i]) * Prb->ntgu[i]
      - (Prb->sig12[i] / Prb->sig11[i]) * Prb->dx[i];
  }

  return 0;
}


/**
  Finds the maximum step size s that can be used in our update
  x_{k+1} = x_k + sdx
  u_{k+1} = u_k + sdu

  that satisfies the constraints that
  f1(x,u) =x - u <=0
  f2(x,u) = -x - u <=0
  f3(x,u) = ||Ax-b||_2^2 - epsilon^2<=0

  In principle, this should be unneccesary, since the line search should reject
  values which don't satisfy the constraints, as they will yield an infinite
  functional. However, this should be faster, that iterating through
  the line search.
 */
double _l1c_l1qc_find_max_step(l1c_l1qcProb *Prb){

  l1c_int m = Prb->m;
  l1c_int n = Prb->n;
  double *Adx = Prb->DWORK7[0];
  double aqe = 0.0, bqe = 0.0, cqe=0.0;
  double smax = 0.0, root = 0.0;
  double min_u1 = 0.0, min_u2=0.0;

  double epsilon = Prb->epsilon;
  l1c_int i=0;
  Prb->ax_funs.Ax(Prb->dx, Adx);  //Adx = A*dx

  aqe = cblas_ddot(n, Adx, 1, Adx, 1);
  bqe = 2.0 * cblas_ddot(n, Prb->r, 1, Adx, 1);
  cqe = cblas_ddot(n, Prb->r, 1, Prb->r, 1) - epsilon * epsilon;

  root = (-bqe + sqrt( bqe*bqe - 4 * aqe * cqe)) / ( 2* aqe);

  /*Towards the end, the root becomes complex. When that happens, cant take
   the real part, because that is just -bqe which is less than zero.*/
  if ( isnan(root) ){
    printf("Warning: maximum step which satisfies cone constraints has become complex.\n");
    printf("Trying smax=1.0 for the cone-constraint portion.\n");
    root = 1.0;
  }

  /*
    We require that the new points  x_{k+1} = x_k + sdx and  u_{k+1} = u_k + sdu
    satisfy:
    f1(x_{k+1},u_{k+1}) = x_k - uk +s(dx -du) <=0
    f2(x_{k+1},u_{k+1}) = -x_k - uk +s(dx -du) <=0

    if dx-du<0, then this will be the case for any positive s, provided
    f1(x_k,u_k) and f2(x_k,u_k) satisfied the constraints. Thus, we only
    care about the cases when,  dx-du>0 or -dx-du>0.
   */
  min_u1 = INFINITY;
  for (i=0; i<m; i++){
    if ( !(Prb->dx[i] - Prb->du[i] > 0) ){
      continue;
    }else{
      min_u1 = min(min_u1, -Prb->fu1[i] / (Prb->dx[i] - Prb->du[i]) );
    }
  }
  min_u2 = INFINITY;
  for (i=0; i<m; i++){
    if ( !(-Prb->dx[i] - Prb->du[i] > 0) ){
      continue;
    }else{
      min_u2 = min(min_u2, -Prb->fu2[i] / (-Prb->dx[i] - Prb->du[i]));
    }
  }

  smax = min(min_u1, min_u2);
  smax = min(smax, root);
  return min(1.0, smax) * 0.99;
}


LSStat _l1c_l1qc_line_search(l1c_int m, double *x, double *u, l1c_l1qcProb *Prb,
                             double g_dot_dxu, LSParams ls_params, double **DWORK5){

  LSStat ls_stat = {.flin=0, .step=0, .status=0};
  int iter=0;
  double flin = 0.0;

  double step = ls_params.s;
  /* Make things accessible */
  double *xp = DWORK5[0];
  double *up = DWORK5[1];

  double f0 = Prb->f_val;

  for (iter=1; iter<=32; iter++){
    /* xp = x + s*dx etc*/
    l1c_daxpy_z(m, step, Prb->dx, x, xp);
    l1c_daxpy_z(m, step, Prb->du, u, up);

    /* Evaluate functional at the new point. */
    _l1c_l1qc_f_eval(Prb, xp, up);

    if ( isnan(Prb->f_val)) {
      step = ls_params.beta * step;
      continue;
    }

    /*
      De-rated linear approximation at the current point. We have
      g_dot_dxu = gradf'*[dx; du]
      So flin = f(x0) + alp * s * [Del x; Del u] * [dx, du]
    */
    flin = f0 + ls_params.alpha * step * g_dot_dxu;

    if (Prb->f_val <= flin){ /* Sufficient decrease test */
      cblas_dcopy(m, xp, 1, x, 1);
      cblas_dcopy(m, up, 1, u, 1);

      ls_stat.flin = flin;
      ls_stat.step = step;
      ls_stat.iter = iter;
      ls_stat.status = 0;
      return ls_stat;
    }

    step = ls_params.beta * step;
  }

  char spc[] = "   ";
  printf("%sBacktracking line search failed, returning previous iterate.\n", spc);
  printf("%sLast line-search values:\n", spc);
  printf("%s                        iter = %d\n", spc, iter);
  printf("%s                        fp   = %.10e\n", spc, Prb->f_val);
  printf("%s                        flin = %.10e\n", spc, flin);
  printf("%s                        fp-flin = %.10f\n", spc, Prb->f_val - flin);
  printf("%s                        step = %.10f\n", spc, step);
  printf("%s                        ls_params.s = %.10f\n", spc, ls_params.s);

  ls_stat.flin = flin;
  ls_stat.step = step;
  ls_stat.iter = min(iter, 32);
  ls_stat.status = 1;
  return ls_stat;
}

int _l1c_l1qc_newton_init(l1c_int m, double *x, double *u,  l1c_L1qcOpts *params){
  l1c_int i = 0;
  double x_max = 0.0, tmp = 0.0;

  /* Initialize u */
  x_max = -INFINITY;
  for (i=0; i<m; i++){
    x_max = max(x_max, fabs(x[i]));
  }
  for (i=0; i<m; i++){
    u[i] = 0.95 * fabs(x[i]) + 0.10 *x_max;
  }

  /* choose initial value of tau so that the duality gap after the first
     step will be about the original norm */
  //tau = max((2*m+1)/sum(abs(x0)), 1);
  tmp = l1c_dnorm1(m, x);
  tmp = (double)(2*m+1) / tmp;
  params->tau = max(tmp, 1);

  // If user has set lbiter, dont compute a different one.
  if (params->lbiter == 0){
    tmp = log (2 * (double)m + 1) - log(params->lbtol) - log(params->tau);
    params->lbiter =(int) ceil (tmp / log(params->mu));
  }
  return 0;
}


int _l1c_l1qc_check_feasible_start(l1c_l1qcProb *Prb, double *x){

  /* Compute Ax - b = r */
  l1c_int n = Prb->n;
  double *b  = Prb->b;
  double epsilon = Prb->epsilon;
  double *DWORK = Prb->DWORK7[0];

  Prb->ax_funs.Ax(x, DWORK);
  cblas_daxpy(Prb->n, -1.0, b, 1, DWORK, 1);

  double tmp = cblas_dnrm2(Prb->n, DWORK, 1)  - epsilon;
  if (tmp > 0){
    // Using minimum-2 norm  x0 = At*inv(AAt)*y.') as they
    // will require updates to cgsolve.
    printf("epsilon = %f,  ||r|| = %.20f\n", epsilon, cblas_dnrm2(n, DWORK, 1));
    return L1C_INFEASIBLE_START;
  }else
    return 0;
}

int _l1c_l1qcProb_new(l1c_l1qcProb *Prb, l1c_int m, l1c_int n, double *b,
                      l1c_L1qcOpts params, l1c_AxFuns ax_funs){
  Prb->m = m;
  Prb->n = n;
  Prb->b = b;
  Prb->epsilon = params.epsilon;
  Prb->tau = params.tau;

  Prb->r = l1c_calloc_double(n);
  Prb->fu1 = l1c_calloc_double(m);
  Prb->fu2 = l1c_calloc_double(m);
  Prb->f_val = 0;
  Prb->fe_val = 0;

  Prb->w1p = l1c_calloc_double(m);
  Prb->dx  = l1c_calloc_double(m);
  Prb->du  = l1c_calloc_double(m);
  Prb->sig11  = l1c_calloc_double(m);
  Prb->sig12  = l1c_calloc_double(m);
  Prb->ntgu   = l1c_calloc_double(m);
  Prb->gradf  = l1c_calloc_double(2*m);
  Prb->DWORK7 = l1c_calloc_double_2D(7, m);

  if ( !Prb->r || !Prb->fu1 || !Prb->fu2
       || !Prb->w1p || !Prb->dx || !Prb->du
       || !Prb->sig11 || !Prb->sig12 ||!Prb->ntgu
       || !Prb->gradf || !Prb->DWORK7){
    return L1C_OUT_OF_MEMORY;
  }

  Prb->ax_funs = ax_funs;
  return 0;

}

void _l1c_l1qcProb_delete(l1c_l1qcProb *Prb){
  /* Do nothing if we got a null pointer. */
  if (Prb){
    l1c_free_double(Prb->r);
    l1c_free_double(Prb->fu1);
    l1c_free_double(Prb->fu2);
    l1c_free_double(Prb->w1p);
    l1c_free_double(Prb->dx);
    l1c_free_double(Prb->du);
    l1c_free_double(Prb->sig11);
    l1c_free_double(Prb->sig12);
    l1c_free_double(Prb->ntgu);
    l1c_free_double(Prb->gradf);
    l1c_free_double_2D(7, Prb->DWORK7);
  }

}



/**
 *
 * @param[in] m size of x.
 * @param[in,out] x result of optimization.
 * @param[in] n size of b
 * @param[in] b The measurements.
 * @param[in] params a struct of options.
 * @param[in] Ax_funs a struct of function pointers.
 */
l1c_LBResult l1c_l1qc_newton(l1c_int m, double *x, l1c_int n, double *b,
                     l1c_L1qcOpts params, l1c_AxFuns Ax_funs){

  // Line search parameters.
  LSStat ls_stat = {.flin=0, .step=0, .status=0};
  LSParams ls_params = { .alpha = 0.01, .beta = 0.5, .s = 1.0};

  l1c_CgParams cg_params = {.tol = params.cg_tol,
                            .max_iter = params.cg_maxiter,
                            .verbose = params.cg_verbose};
  l1c_CgResults cg_results;
  l1c_LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  int iter=0, lb_iter = 0;

  double lambda2 = 0.0, g_dot_dxu=0.0;
  double l1_prev = 0, l1 = 0, rate=0;
  l1c_l1qcProb l1qc_prob = {.m=0, .n=0, .b=NULL, .tau=0, .epsilon=0,
                            .r=NULL, .fu1=NULL, .fu2=NULL,
                            .f_val=0, .fe_val=0,
                            .w1p=NULL, .dx=NULL, .du=NULL, .sig11=NULL,
                            .sig12=NULL, .ntgu=NULL, .gradf=NULL, .DWORK7=NULL};

  double *u   = l1c_calloc_double(m);

  if (!u){
    lb_res.status = L1C_OUT_OF_MEMORY;
    goto exit;
  }
  if( _l1c_l1qcProb_new(&l1qc_prob, m, n, b, params, Ax_funs) ){
    lb_res.status = L1C_OUT_OF_MEMORY;
      goto exit;
  }

  _l1c_l1qc_newton_init(m, x, u, &params);
  l1qc_prob.tau = params.tau;


  if (_l1c_l1qc_check_feasible_start(&l1qc_prob, x) ){
    printf("Starting point is infeasible, exiting\n");
    lb_res.status = L1C_INFEASIBLE_START;
    goto exit;
  }


  if (params.verbose >0) {
    printf("Total Log-Barrier iterations:  %d \n", (int)params.lbiter);
    printf("Original l1-norm: %f, original functional %f\n", l1c_dnorm1(m, x), l1c_dsum(m, u));
  }

  /* ---------------- MAIN **TAU** ITERATION --------------------- */
  for (lb_iter=1; lb_iter<=params.lbiter; lb_iter++){

    /*Re-initialize dx to zero every lb-iteration. */
    l1c_init_vec(m, l1qc_prob.dx, 0.0);

    _l1c_l1qc_f_eval(&l1qc_prob, x, u);

    l1_prev = INFINITY;
    /* ---------------- MAIN Newton ITERATION --------------------- */
    for (iter=1; iter<+params.newton_max_iter; iter++){

      /* Setup warm start for CG solver.
         warm_start_cg ==1 means we use dx itself, ie., a NOP.
       */
      if (params.warm_start_cg==0){
        cblas_dscal(m, 0, l1qc_prob.dx, 1);
      }else if (params.warm_start_cg==2){
        cblas_dscal(m, ls_stat.step, l1qc_prob.dx, 1);
      }
     /* Note: we do not need to re-evaluate the functional because we now do
        that inside the line search. Thus, f, fe, fu1, fu2 and r provided by
        linesearch should be exact given the current stepsize.
      */


      if(_l1c_l1qc_descent_dir(&l1qc_prob, cg_params, &cg_results)){
        printf("Unable to solve system for Newton descent direction.\n");
        printf("Returning previous iterate\n");
      }

      lb_res.total_cg_iter +=cg_results.cgiter;


      /* -------------- Line Search --------------------------- */
      ls_params.s = _l1c_l1qc_find_max_step(&l1qc_prob);

      g_dot_dxu= cblas_ddot(m, l1qc_prob.gradf, 1, l1qc_prob.dx, 1)
        + cblas_ddot(m, l1qc_prob.gradf+m, 1, l1qc_prob.du, 1);

      ls_stat = _l1c_l1qc_line_search(m, x, u, &l1qc_prob, g_dot_dxu, ls_params, l1qc_prob.DWORK7);


      if (ls_stat.status > 0){
        /* When the line search fails, x and u should not get updated, but l1qc_prob was
         evaluated at the last test point, x+s*dx.
         Since we end the inner (newton) iterations, this is fine because at the start of the
         next log-barrier iteration, we evaluate the functional explicetly again.
         */

        printf("Line search failed at newton iteration %d\n", iter);
        break;
      }

      /* ------------------Report and Check for exit --------------*/
      lambda2 = -g_dot_dxu;

      if (params.verbose >1){
        newton_report(iter, m, &l1qc_prob, lambda2, cg_results, ls_stat);
      }

      /*Check for early exit */
      if (lambda2/2 < params.newton_tol){
        break;
      }

      l1 = l1c_dnorm1(m, x);
      rate = fabs(l1_prev - l1)/l1;

      l1_prev = l1;
      if (rate < params.l1_tol){
        if(params.verbose > 0){
          printf("rate = %.9f < %.9f, stopping newton iteration\n", rate, params.l1_tol);
        }
        break;
      }


    }/* Newton iter*/


    /* ----- Update tau ------ */
    lb_res.total_newton_iter += iter;

    params.tau = params.tau * params.mu;
    l1qc_prob.tau = params.tau;

    /* Report on log-barrier progress. */
    if (params.verbose >0){
      lb_report(lb_iter, m, u, l1, params, lb_res);
    }

  } /* LB-iter */

  /* ----- Cleanup -------------------- */


 exit:
  l1c_free_double(u);
  _l1c_l1qcProb_delete(&l1qc_prob);

  lb_res.l1 = l1c_dnorm1(m, x);
  return lb_res;

}/* MAIN ENDS HERE*/



/*
Prints status report for the outer, log-barrier iterations.
 */
static void
lb_report(int lb_iter, int m, double *u, double l1,
          l1c_L1qcOpts params, l1c_LBResult lb_res){

  printf("\n*********************************************************************");
  printf("***********************************\n");
  printf("* LB iter: %d, l1: %.3f, fctl: %8.3e, tau: %8.3e, total newton-iter: %d, Total CG iter=%d *\n",
         lb_iter, l1, l1c_dsum(m, u), params.tau, lb_res.total_newton_iter, lb_res.total_cg_iter);
  printf("***********************************************************************");
  printf("*********************************\n\n");

}


/*
  Prints status reports for the inner Newton iterations.
 */
static void
newton_report(int iter, int m, l1c_l1qcProb *Prb, double lambda2,
              l1c_CgResults cg_results, LSStat ls_status){

  /* want norm [dx; du] */
  double  stepsize = cblas_ddot(m, Prb->dx, 1, Prb->dx, 1) + cblas_ddot(m, Prb->du, 1, Prb->du, 1);
  stepsize = stepsize * ls_status.step;
  if (iter == 1){
    printf("Newton-iter | Functional | Newton decrement |  Stepsize  |  cg-res");
    printf(" | cg-iter | backiter |   s    | \n");
    printf("-----------------------------------------------------------------------");
    printf("----------------------------\n");
  }
  /*            NI         fcnl         dec         sz     cgr       cgI        BI     s  */
  printf("     %3d       %8.3e       %08.3e    % 8.3e   %08.3e     %3d       %2d       %.3g \n",
         (int)iter, Prb->f_val, lambda2/2, stepsize, cg_results.cgres, (int)cg_results.cgiter,
         (int)ls_status.iter, ls_status.step);
}
