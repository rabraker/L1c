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
/* Defined in test_l1c.c*/
extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;


int load_small_data(double **A, double **x, double **b, l1c_int *N, l1c_int *na,
                    l1c_int *max_iter, double *tol){
  l1c_int Nx=0, Nb= 0, status=0;
  char *fpath = fullfile(test_data_dir, "cgsolve_small01.json");
  if(load_file_to_json(fpath, &test_data_json) ){
    ck_abort();
  }
  free(fpath);

  status += extract_json_int(test_data_json, "max_iter", max_iter);
  status += extract_json_double(test_data_json, "tol", tol);

  status +=extract_json_double_array(test_data_json, "x", x, &Nx);
  status +=extract_json_double_array(test_data_json, "b", b, &Nb);
  status += extract_json_double_array(test_data_json, "A", A, na);

  if (status){
    fprintf(stderr, "Error loading json data\n");
    goto fail;
  }
  *N = Nx;

  if ( (Nx != Nb) || ( (l1c_int)(Nb* (Nb +1)/2) != *na) ){
    fprintf(stderr, "Error: Array size mismatch. Aborting\n");
    goto fail;
  }

  return 0;

 fail:
  free_double(*A);
  free_double(*b);
  free_double(*x);
  return 1;

}


START_TEST(test_cgsolve)
{
  double tol =0.0; //= 1e-6;
  l1c_int max_iter;
  CgParams cgp;
  CgResults cgr;

  double *A=NULL, *x=NULL, *x_exp=NULL, *b=NULL, *Dwork=NULL;
  int status = 0;
  l1c_int N=0, na= 0;
  if (load_small_data(&A, &x_exp, &b, &N, &na, &max_iter, &tol)){
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = malloc_double(N);
  Dwork = malloc_double(N*4);
  if (!Dwork || !x){
    status +=1;
    goto exit;
  }

  /* Must initialize x now, with warm starting. */
  for(int i=0; i<N; i++){
    x[i] = 0.0;
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
  char *fpath = fullfile(test_data_dir, "descent_data.json");

  Hess_data h11p_data;
  double *atr=NULL, *sigx=NULL, *dx0=NULL, *dx_exp=NULL;
  double *w1p=NULL, *DWORK_4N=NULL;
  double  fe,cgtol,tau = 0;
  CgResults cgr;
  CgParams cgp = {.verbose=0, .max_iter=0, .tol=0};

  l1c_int N, M, cg_maxiter, status=0;
  l1c_int *pix_idx;

  if (load_file_to_json(fpath, &test_data_json)){
    ck_abort_msg("Error loading data in test_cgsolve_h11p\n");
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
  h11p_data.AtAx = dct_MtEt_EMx;

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

  free(fpath);

  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

  if (status){
    ck_abort();
  }

}
END_TEST


/* Test the matrix multiplication function, Ax_sym.*/
START_TEST(test_cgsolve_Ax_sym){

  char *fpath = fullfile(test_data_dir, "ax_sym.json");

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
    fprintf(stderr, "Error allocating memory\n");
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

  free(fpath);
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
  if (status){
    ck_abort();
  }

}
END_TEST

/* ----------------------------- cgsolve_diag_precond--------------------------------*/

START_TEST(test_cgsolve_diag_precond)
{
  double tol =0.0; //= 1e-6;
  l1c_int max_iter;
  CgParams cgp;
  CgResults cgr;

  double *A=NULL, *x=NULL, *x_exp=NULL, *b=NULL;
  double *Dwork=NULL, *M_inv_diag=NULL;
  int status = 0;
  l1c_int N=0, na= 0;
  if (load_small_data(&A, &x_exp, &b, &N, &na, &max_iter, &tol)){
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = malloc_double(N);
  M_inv_diag = malloc_double(N);
  Dwork = malloc_double(N*5);
  if (!Dwork || !x || !M_inv_diag){
    status +=1;
    goto exit;
  }

  /* Must initialize x now, with warm starting. */
  for(int i=0; i<N; i++){
    x[i] = 0.0;
    M_inv_diag[i] = 1.0;
  }
  cgsolve_diag_precond(x, b, M_inv_diag, N, Dwork, Ax_sym, A, &cgr, cgp);

  ck_assert_double_array_eq_tol(N, x_exp, x, TOL_DOUBLE);
  goto exit;

 exit:
  free_double(A);
  free_double(x);
  free_double(x_exp);
  free_double(b);
  free_double(Dwork);
  free_double(M_inv_diag);

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

  TCase *tc_cgsolve, *tc_cgsolve_precond;
  s = suite_create("cgsolve");

  tc_cgsolve = tcase_create("cgsolve");
  tcase_add_test(tc_cgsolve, test_cgsolve);
  tcase_add_test(tc_cgsolve, test_cgsolve_h11p);
  tcase_add_test(tc_cgsolve, test_cgsolve_Ax_sym);
  suite_add_tcase(s, tc_cgsolve);

  tc_cgsolve_precond = tcase_create("cgsolve_precond");
  tcase_add_test(tc_cgsolve_precond, test_cgsolve_diag_precond);

  suite_add_tcase(s, tc_cgsolve_precond);

  return s;

}
