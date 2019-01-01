/*
This is a test suite for the fourier integration functions contained in
../src/ss_fourier_functions.c. The data used to test these functions is
contained in the header file test_data_ss_ff.h, which defines several global
variables. The header file is generated from the matlab script called
generate_test_data.m

This test suite uses the libcheck framework. On my computer, this got installed
into /usr/local/lib, which was not by default found by the system. Thus, I have
to do
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib


  https://libcheck.github.io/check/doc/check_html/check_4.html#No-Fork-Mode

 */


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <fftw3.h>

#include "l1qc_common.h"

#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"
#include "dct.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"



static cJSON *test_data_json;


START_TEST (test_l1qc_newton_1iter)
{
  CgParams cgp = {.verbose=0, .max_iter=0.0, .tol=0.0};
  NewtParams params;
  char fpath[] = "test_data/lb_test_data_iter_2.json";
  double *x0=NULL, *u=NULL, *b=NULL;
  double *x1_exp=NULL, *u1_exp=NULL;

  LBResult lb_res;
  double epsilon = 0., tau_exp=0., lbtol=0.;
  double mu = 0.0;
  int N=0, M=0, status=0, lbiter_exp=0;
  int *pix_idx=NULL;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_l1qc_newton_1iter\n");
    ck_abort();
  }
  status +=extract_json_double_array(test_data_json, "x0", &x0, &N);
  status +=extract_json_double_array(test_data_json, "xp", &x1_exp, &N);
  status +=extract_json_double_array(test_data_json, "up", &u1_exp, &N);
  status +=extract_json_double_array(test_data_json, "b", &b, &M);

  status +=extract_json_double(test_data_json, "epsilon", &epsilon);
  status +=extract_json_double(test_data_json, "mu", &mu);
  status +=extract_json_double(test_data_json, "lbtol", &lbtol);
  status +=extract_json_double(test_data_json, "newtontol", &params.newton_tol);
  status +=extract_json_int(test_data_json, "newtonmaxiter", &params.newton_max_iter);
  status +=extract_json_double(test_data_json, "cgtol", &cgp.tol);
  status +=extract_json_int(test_data_json, "cgmaxiter", &cgp.max_iter);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);

  status +=extract_json_double(test_data_json, "tau", &tau_exp);

  status +=extract_json_int(test_data_json, "lbiter", &lbiter_exp);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);

  u = malloc_double(N);
  if (status | !u){
    status += 1;
    goto exit1;
  }

  params.verbose = 1;
  params.mu = mu;
  params.lbtol = lbtol;
  params.epsilon = epsilon;
  params.lbiter = 2;
  params.cg_params = cgp;

  lb_res = l1qc_newton(N, x0, u, b, M, pix_idx,params);

  ck_assert_double_array_eq_tol(N, x1_exp, x0,  0.00001);
  ck_assert_double_array_eq_tol(N, u1_exp, u,   0.00001);

  ck_assert_int_eq(0, lb_res.status);


  goto exit1;

 exit1:
  free_double(x0);
  free_double(x1_exp);
  free_double(u1_exp);
  free_double(u);
  free_double(b);
  free(pix_idx);

  dct_destroy();
  free_json_text_data();
  cJSON_Delete(test_data_json);

#ifdef _USEMKL_
  mkl_free_buffers();
#endif


  if (status){
    perror("Error Loading json data in 'test_l1qc_newton_1ter()'. Aborting\n");
    ck_abort();
  }

}
END_TEST

/*------------------------------------------ */
START_TEST (test_newton_init_regres1)
{
  NewtParams params;
  int N=4, M=2;
  int pix_idx[] = {1,2};

  double x[] = {1.0, 2.0, 3.0, 4.0};
  double u[] = {0,0,0,0};
  double u_exp[] = {    1.3500,
                      2.3000,
                      3.2500,
                      4.2000};

  params.lbtol = 0.1;
  params.mu = 0.1;

  newton_init(N, x, u, &params,  M, pix_idx);

  ck_assert_double_array_eq_tol(N, u_exp, u,  TOL_DOUBLE_SUPER);

  dct_destroy();
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST



START_TEST (test_newton_init)
{
  NewtParams params;
  char fpath[] = "test_data/newton_init_data.json";
  double *x=NULL, *u=NULL, *u_exp=NULL, *b=NULL, *Dwork=NULL;
  double epsilon = 0., tau_exp=0., lbtol=0.;
  double mu = 0.0;
  int N=0, M=0, status=0, ret=0, lbiter_exp=0;
  int *pix_idx=NULL;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_newton_init\n");
    ck_abort();
  }
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u_exp, &N);
  status +=extract_json_double_array(test_data_json, "b", &b, &M);

  status +=extract_json_double(test_data_json, "epsilon", &epsilon);
  status +=extract_json_double(test_data_json, "mu", &mu);
  status +=extract_json_double(test_data_json, "lbtol", &lbtol);
  status +=extract_json_double(test_data_json, "tau", &tau_exp);

  status +=extract_json_int(test_data_json, "lbiter", &lbiter_exp);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);
  u = malloc_double(N);
    Dwork = malloc_double(N);
  if (status | !u | !Dwork){
    perror("Error Loading json data in 'test_newton_init()'. Aborting\n");
    goto exit1;
  }

  params.mu = mu;
  params.lbtol = lbtol;
  params.epsilon = epsilon;
  params.lbiter = 0;
  ret= newton_init(N, x, u, &params, M, pix_idx);

  ck_assert_double_array_eq_tol(N, u_exp, u,  TOL_DOUBLE_SUPER);
  ck_assert_double_eq_tol(tau_exp, params.tau,  TOL_DOUBLE_SUPER*100);
  ck_assert_int_eq(lbiter_exp, params.lbiter);
  ck_assert_int_eq(0, ret);

  params.lbiter = 1;
  ret= newton_init(N, x, u, &params,  M, pix_idx);
  ck_assert_int_eq(1, params.lbiter);

  goto exit1;

 exit1:
  free_double(x);
  free_double(u);
  free_double(u_exp);
  free_double(b);
  free(pix_idx);
  free_double(Dwork);

  free_json_text_data();
  cJSON_Delete(test_data_json);

  if (! status){
    goto exit2;
  }else{
#ifdef _USEMKL_
    mkl_free_buffers();
#endif
    ck_abort();
  }
 exit2:
  dct_destroy();
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST

START_TEST (test_find_max_step)
{
  GradData gd;
  char fpath[] = "test_data/find_max_step_data.json";
  double *dx, *du, *Adx, *fu1, *fu2, *r;
  double epsilon, smax, smax_exp = 0.0;
  int N, M, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_find_max_stept\n");
    ck_abort();
  }
  status +=extract_json_double_array(test_data_json, "dx", &dx, &N);
  status +=extract_json_double_array(test_data_json, "du", &du, &N);
  status +=extract_json_double_array(test_data_json, "Adx", &Adx, &N);
  status +=extract_json_double_array(test_data_json, "fu1", &fu1, &N);
  status +=extract_json_double_array(test_data_json, "fu2", &fu2, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &M);
  status +=extract_json_double(test_data_json, "smax", &smax_exp);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);

  if (status){
    perror("Error Loading json data in 'test_find_max_step()'. Aborting\n");
    ck_abort();
  }
  //s = find_max_step(dx, du, Adx, fu1, fu2, r, epsilon, save_on, jopts);
  gd.dx = dx;
  gd.du = du;
  gd.Adx = Adx;

  smax = find_max_step(N, gd, fu1, fu2, M, r, epsilon);
  ck_assert_double_eq_tol(smax_exp, smax,  TOL_DOUBLE*100);

  free_double(dx);
  free_double(du);
  free_double(Adx);
  free_double(fu1);
  free_double(fu2);
  free_double(r);

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_compute_descent)
{
  char fpath[] = "test_data/descent_data.json";

  GradData gd;
  CgParams cgp;
  CgResults cgr;

  double *fu1, *fu2, *atr, *dx_exp, *du_exp;
  double *sig11_exp, *sig12_exp, *w1p_exp;

  double *DWORK_6N;
  double  fe,tau,cgres_exp,cgtol = 0;

  int cg_maxiter, cgiter_exp, N, Npix, status=0;
  int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_compute_descent\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "fu1", &fu1, &N);
  status +=extract_json_double_array(test_data_json, "fu2", &fu2, &N);
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);

  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  status +=extract_json_double(test_data_json, "fe", &fe);
  status +=extract_json_double(test_data_json, "cgtol", &cgtol);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_int(test_data_json, "cgmaxiter", &cg_maxiter);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "dx", &dx_exp, &N);
  status +=extract_json_double_array(test_data_json, "du", &du_exp, &N);
  status +=extract_json_double_array(test_data_json, "sig11", &sig11_exp, &N);
  status +=extract_json_double_array(test_data_json, "sig12", &sig12_exp, &N);
  status +=extract_json_double_array(test_data_json, "w1p", &w1p_exp, &N);

  status +=extract_json_double(test_data_json, "cgres", &cgres_exp);
  status +=extract_json_int(test_data_json, "cgiter", &cgiter_exp);
  if (status){
    perror("Error Loading json data in 'test_compute_descent()'. Aborting\n");
    ck_abort();
  }

  DWORK_6N = malloc_double(6*N);
  gd.w1p = malloc_double(N);
  gd.dx = malloc_double(N);
  gd.du = malloc_double(N);
  gd.gradf = malloc_double(2*N);
  gd.Adx = malloc_double(2*N);
  gd.sig11 = malloc_double(N);
  gd.sig12 = malloc_double(N);
  gd.ntgu = malloc_double(N);
  if ( (!DWORK_6N) | (!gd.w1p) | (!gd.dx) | (!gd.du) | (!gd.gradf) | (!gd.Adx)
       |(!gd.sig11)| (!gd.sig12) |(!gd.ntgu) ){
    perror("Error allocating memory\n");
  }


  cgr.cgres = 0.0;
  cgr.cgiter = 0;
  cgp.verbose = 0;
  cgp.max_iter = cg_maxiter;
  cgp.tol = cgtol;

  /* Setup the DCT */
  dct_setup(N, Npix, pix_idx);

  compute_descent(N, fu1, fu2, atr, fe, tau, gd, DWORK_6N, cgp, &cgr);
  //get_gradient(N, fu1, fu2, DWORK_5N, atr, fe, tau, gd);
  /* ----- Now check: if this fails, look at *relative* tolerance. Here, we check
     absolute tolerance, but for real data,  y_exp is huge 1e22. Should probably normalize
     inputs to 1??
  */

  /*The next three should already be checked by test_get_gradient, but we can
   do it here to, to make sure things are staying sane.*/
  ck_assert_double_array_eq_tol(N, sig11_exp, gd.sig11, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, sig12_exp, gd.sig12, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, w1p_exp, gd.w1p, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, dx_exp, gd.dx, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, du_exp, gd.du, TOL_DOUBLE*100);

  ck_assert_int_eq(cgr.cgiter, cgiter_exp);
  ck_assert_double_eq_tol(cgr.cgres, cgres_exp, TOL_DOUBLE*100);

  /* /\* ----------------- Cleanup --------------- *\/ */
  free_double(fu1);
  free_double(fu2);
  free_double(atr);
  free_double(dx_exp);
  free_double(du_exp);
  free_double(sig11_exp);
  free_double(sig12_exp);
  free_double(w1p_exp);
  free_double(DWORK_6N);
  free_double(gd.w1p);
  free_double(gd.dx);
  free_double(gd.du);
  free_double(gd.gradf);
  free_double(gd.Adx);
  free_double(gd.sig11);
  free_double(gd.sig12);
  free_double(gd.ntgu);
  free(pix_idx);

  dct_destroy();

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST





START_TEST(test_H11pfun)
{
  char fpath[] = "test_data/hp11_fun_data.json";

  Hess_data h11p_data;
  double *atr, *sigx, *z, *y_exp, *y;
  double  fe = 0;

  int N, M, status=0;
  int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_H11pfun\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);
  status +=extract_json_double_array(test_data_json, "sigx", &sigx, &N);
  status +=extract_json_double_array_fftw(test_data_json, "z", &z, &N);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);

  status +=extract_json_double(test_data_json, "fe", &fe);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "y_exp", &y_exp, &N);

  if (status){
    perror("Error Loading json data in 'test_H11pfun()'. Aborting\n");
    ck_abort();
  }
  y = malloc_double(N);
  if (!y){
    perror("Unable to allocate memory\n");
  }

  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = 1.0/(fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;
  h11p_data.Dwork_1N = malloc_double(N);

  /* Setup the DCT */
  dct_setup(N, M, pix_idx);

  H11pfun(N, z, y, &h11p_data);

  /* ----- Now check: if this fails, look at *relative* tolerance. Here, we check
     absolute tolerance, but for real data,  y_exp is huge 1e22. Should probably normalize
     inputs to 1??
  */
  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE*100);

  /* ----------------- Cleanup --------------- */
  free_double(atr);
  free_double(sigx);
  fftw_free(z);
  free(pix_idx);

  free_double(y_exp);
  free_double(y);
  free_double(h11p_data.Dwork_1N);
  dct_destroy();

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST


START_TEST(test_get_gradient)
{
  GradData gd;
  double *fu1, *fu2, *sigx, *sigx_exp, *atr;
  double  *w1p_exp, *sig11_exp, *sig12_exp, *gradf_exp, *ntgu_exp;
  double tau; //loaded

  double fe;
  int N,N2, M, status=0;
  char fpath[] = "test_data/descent_data.json";

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "fu1", &fu1, &N);
  status +=extract_json_double_array(test_data_json, "fu2", &fu2, &N);
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);
  // Expected outputs
  status +=extract_json_double_array(test_data_json, "sigx", &sigx_exp, &N);
  status +=extract_json_double_array(test_data_json, "sig11", &sig11_exp, &N);
  status +=extract_json_double_array(test_data_json, "sig12", &sig12_exp, &M);
  status +=extract_json_double_array(test_data_json, "gradf", &gradf_exp, &N2);
  status +=extract_json_double_array(test_data_json, "ntgu", &ntgu_exp, &N);
  status +=extract_json_double_array(test_data_json, "w1p", &w1p_exp, &N);

  status +=extract_json_double(test_data_json, "fe", &fe);
  status +=extract_json_double(test_data_json, "tau", &tau);

  if (status){
    perror("Error Loading json data in 'test_get_gradient()'. Aborting\n");
    ck_abort();
  }

  gd.sig11 = malloc_double(N);
  gd.sig12 = malloc_double(N);
  gd.gradf = malloc_double(2*N);
  gd.w1p = malloc_double(N);
  gd.ntgu = malloc_double(N);
  sigx = malloc_double(N);

  /*-------------------------------------------- */
    get_gradient(N, fu1, fu2, sigx, atr, fe, tau, gd);


  ck_assert_int_eq(2*N, N2);
  /* ----- Now check -------*/
  ck_assert_double_array_eq_tol(N, sig11_exp, gd.sig11,TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(N, sig12_exp, gd.sig12,TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(N, w1p_exp, gd.w1p,TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(N, sigx_exp, sigx,TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(N, ntgu_exp, gd.ntgu,TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(N2, gradf_exp, gd.gradf,TOL_DOUBLE_SUPER*100);

  /* ----------------- Cleanup --------------- */

  free_double(fu1);
  free_double(fu2);
  free_double(sigx);
  free_double(atr);

  free_double(sigx_exp);
  free_double(sig11_exp);
  free_double(sig12_exp);
  free_double(gradf_exp);
  free_double(ntgu_exp);
  free_double(w1p_exp);

  free_double(gd.sig11);
  free_double(gd.sig12);
  free_double(gd.gradf);
  free_double(gd.ntgu);
  free_double(gd.w1p);

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST



START_TEST(test_line_search)
{
  LSStat ls_stat = {.flx=0, .flu = 0, .flin=0, .step=0, .status=0};
  LSParams ls_params;
  GradData gd;
  double *x, *u, *r, *b, *dx, *du, *gradf;
  double tau, epsilon, alpha, beta, s_init; //loaded

  double *fu1p, *fu2p, fe, f, *DWORK_5N;
  double *xp_exp, *up_exp, *rp_exp;
  double *fu1p_exp, *fu2p_exp, fep_exp, fp_exp;
  double flx_exp=0, flu_exp=0, flin_exp=0;

  //double sm;
  int N,N2, M, status=0;
  int *pix_idx=NULL;
  char fpath[] = "test_data/line_search_data.json";

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_line_search\n");
    ck_abort();
  }

  // Inputs to line_search
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &M);
  status +=extract_json_double_array(test_data_json, "b", &b, &M);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);

  status +=extract_json_double_array(test_data_json, "dx", &dx, &N);
  status +=extract_json_double_array(test_data_json, "du", &du, &N);
  //status +=extract_json_double_array(test_data_json, "Adx", &Adx, &M);
  status +=extract_json_double_array(test_data_json, "gradf", &gradf, &N2);

  status +=extract_json_double(test_data_json, "f", &f);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);
  status +=extract_json_double(test_data_json, "alpha", &alpha);
  status +=extract_json_double(test_data_json, "beta", &beta);
  status +=extract_json_double(test_data_json, "s_init", &s_init);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "xp", &xp_exp, &N);
  status +=extract_json_double_array(test_data_json, "up", &up_exp, &N);
  status +=extract_json_double_array(test_data_json, "rp", &rp_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu1p", &fu1p_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu2p", &fu2p_exp, &N);
  status +=extract_json_double(test_data_json, "fep", &fep_exp);
  status +=extract_json_double(test_data_json, "fp", &fp_exp);

  status +=extract_json_double(test_data_json, "flin", &flin_exp);
  status +=extract_json_double(test_data_json, "flx", &flx_exp);
  status +=extract_json_double(test_data_json, "flu", &flu_exp);

  if (status){
    perror("Error Loading json data in 'test_line_search()'. Aborting\n");
    ck_abort();
  }

  ls_params = (LSParams){.alpha = alpha, .beta = beta,
               .tau = tau, .s = s_init, .epsilon = epsilon};

  gd.dx = dx;
  gd.du = du;
  gd.gradf = gradf;
  // gd.Adx = Adx;

  DWORK_5N = malloc_double(5*N);
  if(!DWORK_5N){
    printf("Allocation failed\n");
  }
  fu1p = malloc_double(N);
  fu2p = malloc_double(N);
  dct_setup(N, M, pix_idx);

  ls_stat = line_search(N, M, x, u, r, b, fu1p, fu2p, gd, ls_params, DWORK_5N, &fe, &f);

  // /* ----- Now check -------*/
  ck_assert_double_eq_tol(fep_exp, fe, TOL_DOUBLE);
  ck_assert_double_eq_tol(fp_exp, f, 1e-7);

  ck_assert_double_eq_tol(flx_exp, ls_stat.flx, TOL_DOUBLE);
  ck_assert_double_eq_tol(flu_exp, ls_stat.flu, TOL_DOUBLE*100);
  ck_assert_double_eq_tol(flin_exp, ls_stat.flin, TOL_DOUBLE);

  ck_assert_double_array_eq_tol(N, xp_exp, x, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, up_exp, u, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu1p_exp, fu1p, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu2p_exp, fu2p, TOL_DOUBLE);

  ck_assert_double_array_eq_tol(M, rp_exp, r, TOL_DOUBLE);


  /* ----------------- Cleanup --------------- */
  free_double(x);
  free_double(u);
  free_double(r);
  free_double(dx);
  free_double(du);
  // free_double(Adx);
  free_double(gradf);

  free_double(xp_exp);
  free_double(up_exp);
  free_double(rp_exp);

  free_double(fu1p_exp);
  free_double(fu2p_exp);
  free_double(fu1p);
  free_double(fu2p);

  free_double(DWORK_5N);

  dct_destroy();

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST




START_TEST(test_f_eval)
{
  double *x, *u, *r, tau, epsilon; //loaded
  double *fu1, *fu2, fe, f;
  double *fu1_exp, *fu2_exp, fe_exp, f_exp;

  int N=0, M=0, status=0;
  char fpath[] = "test_data/f_eval_data.json";
  if (  load_file_to_json(fpath, &test_data_json) ){
    perror("Error loading data for test_f_eval\n");
    ck_abort();
  }
  // Inputs to f_eval
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &M);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "fu1_exp", &fu1_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu2_exp", &fu2_exp, &N);
  status +=extract_json_double(test_data_json, "f_exp", &f_exp);
  status +=extract_json_double(test_data_json, "fe_exp", &fe_exp);
  if (status){
    perror("Error Loading json data in 'test_f_eval'. Aborting\n");
    ck_abort();
  }

  fu1 = malloc_double(N);
  fu2 = malloc_double(N);

  f_eval(N,  x, u, M, r, tau, epsilon, fu1, fu2, &fe, &f);

  /* ----- Now check -------*/
  ck_assert_double_eq_tol(fe_exp, fe, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu1_exp, fu1, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu2_exp, fu2, TOL_DOUBLE);

  ck_assert_double_eq_tol(f_exp, f, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */
  free_double(x);
  free_double(u);
  free_double(r);

  free_double(fu1_exp);
  free_double(fu2_exp);

  free_double(fu1);
  free_double(fu2);

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST

START_TEST(test_sum_vec)
{
  int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double sum_exp = 21.0;
  double sum_x = sum_vec(N, x);

  ck_assert_double_eq_tol(sum_exp, sum_x, TOL_DOUBLE);
}
END_TEST


START_TEST(test_logsum)
{
  int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double alpha = 3.0;
  double logsum_x_exp = 13.170924944018758;
  double logsum_x = logsum(N, alpha, x);

  ck_assert_double_eq_tol(logsum_x_exp, logsum_x, TOL_DOUBLE);
  // setup_vectors_SC();
}
END_TEST




Suite *l1qc_newton_suite(void)
{
  Suite *s;

  TCase *tc_linesearch, *tc_newton_init, *tc_feval, *tc_mathfuns, *tc_gradient;
  TCase *tc_l1qc_newton;

  s = suite_create("l1qc_newton");

  tc_l1qc_newton = tcase_create("l1qc_newton");
  tcase_add_test(tc_l1qc_newton, test_l1qc_newton_1iter);

  tc_newton_init = tcase_create("newton_init");
  tcase_add_test(tc_newton_init, test_newton_init);
  tcase_add_test(tc_newton_init, test_newton_init_regres1);

  tc_feval = tcase_create("f_eval");
  tcase_add_test(tc_feval, test_f_eval);

  tc_linesearch = tcase_create("l1qc_linesearch");
  tcase_add_test(tc_linesearch, test_find_max_step);
  tcase_add_test(tc_linesearch, test_line_search);

  tc_gradient = tcase_create("l1qc_gradient");
  tcase_add_test(tc_gradient, test_H11pfun);
  tcase_add_test(tc_gradient, test_get_gradient);
  tcase_add_test(tc_gradient, test_compute_descent);

  tc_mathfuns = tcase_create("l1qc_math_funs");
  tcase_add_test(tc_mathfuns, test_sum_vec);
  tcase_add_test(tc_mathfuns, test_logsum);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_l1qc_newton);
  suite_add_tcase(s, tc_newton_init);
  suite_add_tcase(s, tc_feval);
  suite_add_tcase(s, tc_linesearch);
  suite_add_tcase(s, tc_gradient);
  suite_add_tcase(s, tc_mathfuns);

  return s;

}
