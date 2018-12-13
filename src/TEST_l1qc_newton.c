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

#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"


/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"
// #include "test_data_ss_ff.h"




static cJSON *test_data_json;


START_TEST(test_get_gradient)
{
  GradData gd;
  double *fu1, *fu2, *sigx, *sigx_exp, *atr;
  double  *w1p_exp, *sig11_exp, *sig12_exp, *gradf_exp, *ntgu_exp;
  double tau; //loaded

  double fe;
  int N,N2, M, status=0;
  char fpath[] = "test_data/descent_data.json";

  load_file_to_json(fpath, &test_data_json);

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
  load_file_to_json(fpath, &test_data_json);
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
  load_file_to_json(fpath, &test_data_json);
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
  // setup_vectors_SC();
  // setup_params_SC();
  int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double sum_exp = 21.0;
  double sum_x = sum_vec(N, x);

  ck_assert_double_eq_tol(sum_exp, sum_x, TOL_DOUBLE);
  // setup_vectors_SC();
}
END_TEST

START_TEST(test_log_vec)
{
  // setup_vectors_SC();
  // setup_params_SC();
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

  tcase_add_test(tc_core, test_get_gradient);
  tcase_add_test(tc_core, test_line_search);
  tcase_add_test(tc_core, test_sum_vec);
  tcase_add_test(tc_core, test_logsum);
  tcase_add_test(tc_core, test_log_vec);
  tcase_add_test(tc_core, test_f_eval);



  suite_add_tcase(s, tc_core);

  return s;

}
