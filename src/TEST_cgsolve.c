/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */

#define CK_FLOATING_DIG 15

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

// #include "test_data.h"
#include "cgsolve.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"
#include "dct.h"


#include "check_utils.h"


int load_small_data(double **A, double **x, double **b, int *N,
                    int *max_iter, double *tol){

  cJSON *test_data_json;

  if(load_file_to_json("test_data/cgsolve_small01.json", &test_data_json) ){
    return 1;
  }

  if( extract_json_int(test_data_json, "max_iter", max_iter) )
    return 1;
  if ( extract_json_double(test_data_json, "tol", tol) )
    return 1;

  int Nx, Nb, Na = 0;
  if (extract_json_double_array(test_data_json, "x", x, &Nx) ){
    perror("Error Loading x\n");
    return 1;
  }

  if (extract_json_double_array(test_data_json, "b", b, &Nb) ){
    perror("Error Loading b\n");
    goto end0;
  }

  if (extract_json_double_array(test_data_json, "A", A, &Na) ){
    perror("Error Loading A\n");
    goto end1;
  }

  *N = Nx;
  printf("Nx = %d\n", *N);
  // print_vec(*N, *x, "x");

  if ( (Nx != Nb) || (Nb*Nb != Na) ){
    perror("Error: Array size mismatch. Aborting\n");
    goto end2; // We allocated all, but their sizes don't match.
  }

  return 0;

 end2:
  free(*A);
  goto end1;
 end1:
  free(*b);
  goto end0;
 end0:
  free(*x);
  return 1;

}


START_TEST(test_cgsolve)
{
  double tol =0.0; //= 1e-6;
  int max_iter;
  CgParams cgp;
  CgResults cgr;

  double *A, *x, *x_exp, *b, *Dwork;
  int N, i = 0;
  if (load_small_data(&A, &x_exp, &b, &N, &max_iter, &tol)){
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = calloc(N, sizeof(double));
  Dwork = calloc(N*4, sizeof(double));

 printf("max-iter = %d, tol=%f\n", max_iter, tol);
  cgsolve(x, b, N, Dwork, Ax, A, &cgr, cgp);

  for (i=0; i<N; i++){
    ck_assert_double_eq_tol(x_exp[i], x[i], TOL_DOUBLE);
  }

  free(A);
  free(x);
  free(x_exp);
  free(b);
  free(Dwork);

}
END_TEST

START_TEST(test_cgsolve_h11p){
  cJSON *test_data_json;

  char fpath[] = "test_data/descent_data.json";

  Hess_data h11p_data;
  double *atr, *sigx, *dx, *dx_exp, *w1p, *DWORK_4N;
  double  fe,cgtol,tau = 0;
  CgResults cgr;
  CgParams cgp;

  int N, M, cg_maxiter, status=0;
  int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);
  status +=extract_json_double_array(test_data_json, "sigx", &sigx, &N);
  status +=extract_json_double_array(test_data_json, "w1p", &w1p, &N);
  status +=extract_json_double_array(test_data_json, "dx", &dx_exp, &N);

  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);

  status +=extract_json_double(test_data_json, "fe", &fe);
  status +=extract_json_double(test_data_json, "cgtol", &cgtol);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_int(test_data_json, "cgmaxiter", &cg_maxiter);

  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = 1.0 / (fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;

  dx = calloc(N, sizeof(double));
  if (!dx){
    perror("error allocating memory\n");
  }
  DWORK_4N = calloc(4*N, sizeof(double));
  if (!DWORK_4N){
    perror("error allocating memory\n");
  }

  dct_setup(N, M, pix_idx);
  cgp.max_iter = cg_maxiter;
  cgp.tol = cgtol;
  cgsolve(dx, w1p, N, DWORK_4N, H11pfun, &h11p_data, &cgr, cgp);
  printf("cgres = %f, cgiter = %d\n", cgr.cgres, cgr.cgiter);

  ck_assert_double_array_eq_tol(N, dx_exp, dx, TOL_DOUBLE*100);

  dct_destroy();
}
END_TEST


/* Add all the test cases to our suite
 */
Suite *cgsolve_suite(void)
{
  Suite *s;

  TCase *tc_core;
  s = suite_create("cgsolve");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_cgsolve);
  tcase_add_test(tc_core, test_cgsolve_h11p);

  suite_add_tcase(s, tc_core);

  return s;

}
