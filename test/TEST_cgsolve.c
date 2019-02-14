/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */
#include "config.h"

#define CK_FLOATING_DIG 20

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

// #include "test_data.h"
#include "cgsolve.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"
#include "l1qc_newton.h"
#include "dct.h"

#include "l1c_common.h"
#include "check_utils.h"

cJSON *test_data_json;

int load_small_data(double **A, double **x, double **b, l1c_int *N, l1c_int *na,
                    l1c_int *max_iter, double *tol){

  if(load_file_to_json("test_data/cgsolve_small01.json", &test_data_json) ){
    return 1;
  }

  if( extract_json_int(test_data_json, "max_iter", max_iter) )
    return 1;
  if ( extract_json_double(test_data_json, "tol", tol) )
    return 1;

  l1c_int Nx=0, Nb= 0;
  if (extract_json_double_array(test_data_json, "x", x, &Nx) ){
    perror("Error Loading x\n");
    return 1;
  }

  if (extract_json_double_array(test_data_json, "b", b, &Nb) ){
    perror("Error Loading b\n");
    goto end0;
  }

  if (extract_json_double_array(test_data_json, "A", A, na) ){
    perror("Error Loading A\n");
    goto end1;
  }

  *N = Nx;

  if ( (Nx != Nb) || ( (l1c_int)(Nb* (Nb +1)/2) != *na) ){
    perror("Error: Array size mismatch. Aborting\n");
    goto end2; // We allocated all, but their sizes don't match.
  }

  return 0;

 end2:
  free_double(*A);
  goto end1;
 end1:
  free_double(*b);
  goto end0;
 end0:
  free_double(*x);
  return 1;

}


START_TEST(test_cgsolve)
{
  double tol =0.0; //= 1e-6;
  l1c_int max_iter;
  CgParams cgp;
  CgResults cgr;

  double *A, *x, *x_exp, *b, *Dwork;
  int status = 0;
  l1c_int N=0, na= 0;
  if (load_small_data(&A, &x_exp, &b, &N, &na, &max_iter, &tol)){
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = malloc_double(N);
  /* Must initialize x now, with warm starting. */
  for(int i=0; i<N; i++){
    x[i] = 0.0;
  }
  Dwork = malloc_double(N*4);
  if (!Dwork){
    status +=1;
    goto exit;
  }
  cgsolve(x, b, N, Dwork, Ax_sym, A, &cgr, cgp);

  ck_assert_double_array_eq_tol(N, x_exp, x, TOL_DOUBLE);
  goto exit;

 exit:
  free_double(A);
  free_double(x);
  free_double(x_exp);
  free_double(b);
  free_double(Dwork);


  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
  if (status){
    ck_abort();
  }

}
END_TEST

START_TEST(test_cgsolve_h11p){

  char fpath[] = "test_data/descent_data.json";

  Hess_data h11p_data;
  double *atr=NULL, *sigx=NULL, *dx0=NULL, *dx_exp=NULL;
  double *w1p=NULL, *DWORK_4N=NULL;
  double  fe,cgtol,tau = 0;
  CgResults cgr;
  CgParams cgp = {.verbose=0, .max_iter=0, .tol=0};

  l1c_int N, M, cg_maxiter, status=0;
  l1c_int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "atr", &atr, &N);
  status +=extract_json_double_array(test_data_json, "sigx", &sigx, &N);
  status +=extract_json_double_array(test_data_json, "w1p", &w1p, &N);
  status +=extract_json_double_array(test_data_json, "dx", &dx_exp, &N);
  status +=extract_json_double_array(test_data_json, "dx0", &dx0, &N);

  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &M);

  status +=extract_json_double(test_data_json, "fe", &fe);
  status +=extract_json_double(test_data_json, "cgtol", &cgtol);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_int(test_data_json, "cgmaxiter", &cg_maxiter);


  h11p_data.one_by_fe = 1.0/fe;
  h11p_data.one_by_fe_sqrd = 1.0 / (fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;
  h11p_data.Dwork_1N = malloc_double(N);
  h11p_data.AtAx = dct_MtEt_EMx_new;

  DWORK_4N = malloc_double(4*N);
  double *dx_by_nrm = malloc_double(4*N);
  double *dx_by_nrm_exp = malloc_double(4*N);
  if ( status || !DWORK_4N || !dx_by_nrm || !dx_by_nrm_exp){
    perror("error allocating memory or reading JSON data in test_cgsolve_h11p\n");
    status +=1;
    goto exit;
  }

  dct_setup(N, M, pix_idx);
  cgp.max_iter = cg_maxiter;
  cgp.tol = cgtol;
  cgp.verbose = 0;

  cgsolve(dx0, w1p, N, DWORK_4N, H11pfun, &h11p_data, &cgr, cgp);
  double nrm_exp = cblas_dnrm2(N, dx_exp, 1);
  double nrm = cblas_dnrm2(N, dx0, 1);
  printf("M=%d, norm exp: %.10f, nrm: %.10f\n", (int)M, nrm_exp, nrm);
  for (int i=0; i<N ; i++){
    dx_by_nrm_exp[i] = dx_exp[i]/nrm_exp;
    dx_by_nrm[i] = dx0[i]/nrm;
  }
  ck_assert_double_array_eq_tol(N, dx_by_nrm_exp, dx_by_nrm, TOL_DOUBLE*10);

  goto exit;

 exit:
  free_double(dx_by_nrm_exp);
  free_double(dx_by_nrm);
  free_double(atr);
  free_double(sigx);
  free_double(w1p);
  free_double(dx_exp);
  free(pix_idx);
  free_double(dx0);
  free_double(DWORK_4N);
  free_double(h11p_data.Dwork_1N);
  dct_destroy();


  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

  if (status){
    ck_abort();
  }

}
END_TEST



START_TEST(test_cgsolve_Ax_sym){

  char fpath[] = "test_data/ax_sym.json";

  double *A, *x, *y_exp, *y;

  l1c_int N=0, na=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status +=extract_json_double_array(test_data_json, "A", &A, &na);
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "y", &y_exp, &N);

  y = malloc_double(N);
  if ( (!y) | status){
    perror("error allocating memory\n");
    status +=1;
    goto exit;
  }

  Ax_sym(N, x, y, A);

  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE_SUPER*10);

  // Now, y should be non-zero. Should get the same answer. Regression against having beta !=0,
  // because dspmv computes alpha*A*x + b*y
  Ax_sym(N, x, y, A);
  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE_SUPER*10);

  goto exit;

 exit:
  free_double(A);
  free_double(x);
  free_double(y_exp);
  free_double(y);


  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
  if (status){
    ck_abort();
  }

}
END_TEST

/* Add all the test cases to our suite
 */
Suite *cgsolve_suite(void)
{
  Suite *s;

  TCase *tc_small, *tc_hp11, *tc_Ax;
  s = suite_create("cgsolve");
  tc_small = tcase_create("cg_small");
  tc_hp11 = tcase_create("cg_hp11");
  tc_Ax = tcase_create("cg_Ax");

  tcase_add_test(tc_small, test_cgsolve);
  tcase_add_test(tc_hp11, test_cgsolve_h11p);
  tcase_add_test(tc_Ax, test_cgsolve_Ax_sym);

  suite_add_tcase(s, tc_small);
  suite_add_tcase(s, tc_hp11);
  suite_add_tcase(s, tc_Ax);

  return s;

}
