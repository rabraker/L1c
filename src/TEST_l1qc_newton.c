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


#define CK_FLOATING_DIG 15
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <fftw3.h>

#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"
#include "dct.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"
// #include "test_data_ss_ff.h"

static cJSON *test_data_json;



START_TEST (test_find_max_step)
{
  GradData gd;
  char fpath[] = "test_data/find_max_step_data.json";
  double *dx, *du, *Adx, *fu1, *fu2, *r;
  double epsilon, smax, smax_exp = 0.0;
  int N, M, status=0;
  int *Iwork_2N;
  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
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
  Iwork_2N = calloc(2*N, sizeof(int));
  smax = find_max_step(N, gd, fu1, fu2, M, r, epsilon, Iwork_2N);

  ck_assert_double_eq_tol(smax_exp, smax,  TOL_DOUBLE*100);

  free(dx);
  free(du);
  free(Adx);
  free(fu1);
  free(fu2);
  free(r);
  free(Iwork_2N);
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

  double *DWORK_5N;
  double  fe,tau,cgres_exp,cgtol = 0;

  int cg_maxiter, cgiter_exp, N, Npix, status=0;
  int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
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

  printf("N = %d, Npix=%d\n", N, Npix);
  DWORK_5N = calloc(5*N, sizeof(double));
  if (!DWORK_5N){
    perror("error allocating memory\n");
  }
  gd.w1p = calloc(N, sizeof(double));
  if (!gd.w1p ){
    perror("error allocating memory\n");
  }
  gd.dx = calloc(N, sizeof(double));
  if (!gd.dx){
    perror("error allocating memory\n");
  }
  gd.du = calloc(N, sizeof(double));
  if (!gd.du ){
    perror("error allocating memory\n");
  }
  gd.gradf = calloc(2*N, sizeof(double));
  if (!gd.gradf){
    perror("error allocating memory\n");
  }
  gd.Adx = calloc(2*N, sizeof(double));
  if (!gd.Adx){
    perror("error allocating memory\n");
  }
  gd.sig11 = calloc(N, sizeof(double));
  if (!gd.sig11){
    perror("error allocating memory\n");
  }
  gd.sig12 = calloc(N, sizeof(double));
  if (!gd.sig12 ){
    perror("error allocating memory\n");
  }
  //  gd.w1p = calloc(N, sizeof(double));
  gd.ntgu = calloc(N, sizeof(double));
  if (!gd.ntgu ){
    perror("error allocating memory\n");
  }


  cgr.cgres = 0.0;
  cgr.cgiter = 0;
  cgp.verbose = 0;
  cgp.max_iter = cg_maxiter;
  cgp.tol = cgtol;

  /* Setup the DCT */
  dct_setup(N, Npix, pix_idx);

  compute_descent(N, fu1, fu2, atr, fe, tau, gd, DWORK_5N, cgp, &cgr);
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

  ck_assert_int_eq(cgr.cgiter+1, cgiter_exp);
  ck_assert_double_eq_tol(cgr.cgres, cgres_exp, TOL_DOUBLE*100);

  /* /\* ----------------- Cleanup --------------- *\/ */
  free(fu1);
  free(fu2);
  free(atr);
  free(dx_exp);
  free(du_exp);
  free(sig11_exp);
  free(sig12_exp);
  free(w1p_exp);
  free(DWORK_5N);
  free(gd.w1p);
  free(gd.dx);
  free(gd.du);
  free(gd.gradf);
  free(gd.Adx);
  free(gd.sig11);
  free(gd.sig12);
  free(gd.ntgu);
  free(pix_idx);

  dct_destroy();
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
    perror("Error loading data in test_get_gradient\n");
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
    perror("Error Loading json data in 'test_get_gradient()'. Aborting\n");
    ck_abort();
  }
  y = calloc(N, sizeof(double));
  if (!y){
    perror("Unable to allocate memory\n");
  }

  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = 1.0/(fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;

  /* Setup the DCT */
  dct_setup(N, M, pix_idx);

  H11pfun(N, z, y, &h11p_data);

  /* ----- Now check: if this fails, look at *relative* tolerance. Here, we check
     absolute tolerance, but for real data,  y_exp is huge 1e22. Should probably normalize
     inputs to 1??
  */
  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE*100);

  /* ----------------- Cleanup --------------- */
  free(atr);
  free(sigx);
  fftw_free(z);
  free(y_exp);
  free(y);

  free(pix_idx);
  dct_destroy();
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
  status +=extract_json_double_array(test_data_json, "sigx", &sigx_exp, &N);
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);
  // Expected outputs
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

  gd.sig11 = calloc(N, sizeof(double));
  gd.sig12 = calloc(N, sizeof(double));
  gd.gradf = calloc(2*N, sizeof(double));
  gd.w1p = calloc(N, sizeof(double));
  gd.ntgu = calloc(N, sizeof(double));
  sigx = calloc(N, sizeof(double));

  printf("N = %d, N2=%d \n", N, N2);

  get_gradient(N, fu1, fu2, sigx, atr, fe, tau, gd);

  ck_assert_int_eq(2*N, N2);
  /* ----- Now check -------*/
  ck_assert_double_array_eq_tol(N, sig11_exp, gd.sig11,TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, sig12_exp, gd.sig12,TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, w1p_exp, gd.w1p,TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, sigx_exp, sigx,TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N, ntgu_exp, gd.ntgu,TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(N2, gradf_exp, gd.gradf,TOL_DOUBLE*100);

  /* ----------------- Cleanup --------------- */

  free(fu1);
  free(fu2);
  free(sigx);
  free(sigx_exp);
  free(atr);

  free(gd.sig11);
  free(sig11_exp);
  free(gd.sig12);
  free(sig12_exp);
  free(gd.gradf);
  free(gradf_exp);
  free(ntgu_exp);
  free(gd.ntgu);
  free(w1p_exp);
  free(gd.w1p);


}
END_TEST



START_TEST(test_line_search)
{
  LSParams ls_params;
  GradData gd;
  double *x, *u, *r, *dx, *du, *Adx, *gradf;
  double tau, epsilon, alpha, beta, s; //loaded

  double *fu1p, *fu2p, fe, f, *DWORK_5N;
  double *xp_exp, *up_exp, *rp_exp;
  double *fu1p_exp, *fu2p_exp, fep_exp, fp_exp;
  //double sm;
  int N,N2, M, status=0;
  char fpath[] = "test_data/line_search_data.json";

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to f_eval
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &M);
  printf("---- M-rp = %d\n", M);
  status +=extract_json_double_array(test_data_json, "dx", &dx, &N);
  status +=extract_json_double_array(test_data_json, "du", &du, &N);
  status +=extract_json_double_array(test_data_json, "Adx", &Adx, &M);
  status +=extract_json_double_array(test_data_json, "gradf", &gradf, &N2);

  status +=extract_json_double(test_data_json, "f", &f);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);
  status +=extract_json_double(test_data_json, "alpha", &alpha);
  status +=extract_json_double(test_data_json, "beta", &beta);
  status +=extract_json_double(test_data_json, "s", &s);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "xp", &xp_exp, &N);
  status +=extract_json_double_array(test_data_json, "up", &up_exp, &N);
  status +=extract_json_double_array(test_data_json, "rp", &rp_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu1p", &fu1p_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu2p", &fu2p_exp, &N);
  status +=extract_json_double(test_data_json, "fep", &fep_exp);
  status +=extract_json_double(test_data_json, "fp", &fp_exp);

  if (status){
    perror("Error Loading json data in 'test_line_search()'. Aborting\n");
    ck_abort();
  }

  ls_params.alpha = alpha;
  ls_params.beta = beta;
  ls_params.tau = tau;
  ls_params.s = s;
  ls_params.epsilon = epsilon;

  gd.dx = dx;
  gd.du = du;
  gd.gradf = gradf;
  gd.Adx = Adx;
  printf("N = %d, N2=%d \n", N, N2);
  DWORK_5N = calloc(5*N, sizeof(double));
  if(!DWORK_5N){
    printf("Allocation failed\n");
  }
  fu1p = calloc(N, sizeof(double));
  fu2p = calloc(N, sizeof(double));

  line_search(N, M, x, u, r, fu1p, fu2p, gd, ls_params, DWORK_5N, &fe, &f);

  // /* ----- Now check -------*/
  ck_assert_double_eq_tol(fep_exp, fe, TOL_DOUBLE);
  ck_assert_double_eq_tol(fp_exp, f, TOL_DOUBLE);

  ck_assert_double_array_eq_tol(N, xp_exp, x, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, up_exp, u, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu1p_exp, fu1p, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu2p_exp, fu2p, TOL_DOUBLE);

  ck_assert_double_array_eq_tol(M, rp_exp, r, TOL_DOUBLE);


  /* ----------------- Cleanup --------------- */
  free(x);
  free(u);
  free(r);
  free(dx);
  free(du);
  free(Adx);
  free(gradf);

  free(xp_exp);
  free(up_exp);
  free(rp_exp);

  free(fu1p_exp);
  free(fu2p_exp);
  free(fu1p);
  free(fu2p);

  free(DWORK_5N);
}
END_TEST




START_TEST(test_f_eval)
{
  double *x, *u, *r, tau, epsilon; //loaded
  double *fu1, *fu2, fe, f, *Dwork;
  double *fu1_exp, *fu2_exp, fe_exp, f_exp;

  int N, status=0;
  char fpath[] = "test_data/f_eval_data.json";
  if (  load_file_to_json(fpath, &test_data_json) ){
    perror("Error loading data for test_f_eval\n");
    ck_abort();
  }
  // Inputs to f_eval
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &N);
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
  Dwork = calloc(4*N, sizeof(double));
  fu1 = calloc(N, sizeof(double));
  fu2 = calloc(N, sizeof(double));

  f_eval(N, r, x, u, tau, epsilon, fu1, fu2, &fe, &f, Dwork);

  /* ----- Now check -------*/
  ck_assert_double_eq_tol(fe_exp, fe, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu1_exp, fu1, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(N, fu2_exp, fu2, TOL_DOUBLE);

  ck_assert_double_eq_tol(f_exp, f, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */
  free(x);
  free(u);
  free(r);

  free(fu1_exp);
  free(fu2_exp);

  free(Dwork);
  free(fu1);
  free(fu2);
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

START_TEST(test_log_vec)
{

  int i, N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double alpha = 3.0;
  double log_alp_x [] = {1.098612288668110,   1.791759469228055,   2.197224577336220,
                         2.484906649788000,   2.708050201102210,   2.890371757896165};
  double y[6];

  log_vec(N, alpha, x, y);

  for (i=0; i<N; i++){
    ck_assert_double_eq_tol(log_alp_x[i], y[i], TOL_DOUBLE);
  }
  // setup_vectors_SC();
}
END_TEST

START_TEST(test_logsum)
{
  int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double alpha = 3.0;
  double logsum_x_exp = 13.170924944018758;
  double logsum_x = logsum(N, x, alpha);

  ck_assert_double_eq_tol(logsum_x_exp, logsum_x, TOL_DOUBLE);
  // setup_vectors_SC();
}
END_TEST




Suite *l1qc_newton_suite(void)
{
  Suite *s;

  TCase *tc_core;
  s = suite_create("l1qc_newton");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_find_max_step);
  tcase_add_test(tc_core, test_compute_descent);
  tcase_add_test(tc_core, test_H11pfun);
  tcase_add_test(tc_core, test_get_gradient);
  tcase_add_test(tc_core, test_line_search);
  tcase_add_test(tc_core, test_sum_vec);
  tcase_add_test(tc_core, test_logsum);
  tcase_add_test(tc_core, test_log_vec);
  tcase_add_test(tc_core, test_f_eval);



  suite_add_tcase(s, tc_core);

  return s;

}
