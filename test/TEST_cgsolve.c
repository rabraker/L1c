/*
  Tests for the conjugate gradient solver.

 */
#include "config.h"

#define CK_FLOATING_DIG 20

#include "cblas.h"
#include <check.h>
#include <cjson/cJSON.h>
#include <math.h> //Constants
#include <stdio.h>
#include <stdlib.h>

#include "l1c.h"

#include "json_utils.h"
#include "l1qc_newton.h"
#include "test_constants.h"

#include "check_utils.h"
#include "l1c_math.h"

static cJSON* test_data_json;

/* Defined in test_l1c.c*/
extern char* fullfile(char* base_path, char* name);
extern char* test_data_dir;

/*
   For test routines.
   Computes the matrix-vector product y = A * b, for a symmetric matrix A.
   This is a wrapper for cblas_dspmv.
*/
static void Ax_sym(l1c_int n, double* x, double* b, void* AX_data) {

  double* A = (double*)AX_data;

  cblas_dspmv(CblasRowMajor, CblasUpper, n, 1.0, A, x, 1, 0.0, b, 1);
}

static int load_small_data(
    double** A, double** x, double** b, l1c_int* N, l1c_int* na, l1c_int* max_iter, double* tol) {
  l1c_int Nx = 0, Nb = 0, status = 0;
  char* fpath = fullfile(test_data_dir, "cgsolve_small01.json");
  if (load_file_to_json(fpath, &test_data_json)) {
    ck_abort();
  }
  free(fpath);

  status += extract_json_int(test_data_json, "max_iter", max_iter);
  status += extract_json_double(test_data_json, "tol", tol);

  status += extract_json_double_array(test_data_json, "x", x, &Nx);
  status += extract_json_double_array(test_data_json, "b", b, &Nb);
  status += extract_json_double_array(test_data_json, "A", A, na);

  if (status) {
    fprintf(stderr, "Error loading json data\n");
    goto fail;
  }
  *N = Nx;

  if ((Nx != Nb) || ((l1c_int)(Nb * (Nb + 1) / 2) != *na)) {
    fprintf(stderr, "Error: Array size mismatch. Aborting\n");
    goto fail;
  }

  return 0;

fail:
  l1c_free_double(*A);
  l1c_free_double(*b);
  l1c_free_double(*x);
  return 1;
}

START_TEST(test_cgsolve) {
  double tol = 0.0; //= 1e-6;
  l1c_int max_iter;
  l1c_CgParams cgp;
  l1c_CgResults cgr;

  double *A = NULL, *x = NULL, *x_exp = NULL, *b = NULL, **Dwork = NULL;
  int status = 0;
  l1c_int N = 0, na = 0;
  if (load_small_data(&A, &x_exp, &b, &N, &na, &max_iter, &tol)) {
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = l1c_malloc_double(N);
  Dwork = l1c_malloc_double_2D(4, N);
  if (!Dwork || !x) {
    status += 1;
    goto exit;
  }

  /* Must initialize x now, with warm starting. */
  for (int i = 0; i < N; i++) {
    x[i] = 0.0;
  }

  l1c_cgsolve(N, x, b, Dwork, Ax_sym, A, &cgr, cgp);

  ck_assert_double_array_eq_tol(N, x_exp, x, TOL_DOUBLE);
  goto exit;

exit:
  l1c_free_double(A);
  l1c_free_double(x);
  l1c_free_double(x_exp);
  l1c_free_double(b);
  l1c_free_double_2D(4, Dwork);

  cJSON_Delete(test_data_json);

  if (status) {
    ck_abort();
  }
}
END_TEST

START_TEST(test_cgsolve_h11p) {
  char* fpath = fullfile(test_data_dir, "l1qc_data.json");

  Hess_data h11p_data;
  double *atr = NULL, *sigx = NULL, *dx0 = NULL, *dx_exp = NULL;
  double *w1p = NULL, **DWORK4 = NULL;
  double fe, cgtol, tau = 0;
  l1c_CgResults cgr;
  l1c_CgParams cgp = {.verbose = 0, .max_iter = 0, .tol = 0};

  l1c_int m, n, cg_maxiter, status = 0;
  l1c_int* pix_idx;
  l1c_AxFuns ax_funs;

  if (load_file_to_json(fpath, &test_data_json)) {
    ck_abort_msg("Error loading data in test_cgsolve_h11p\n");
  }

  // Inputs to get_gradient
  status += extract_json_double_array(test_data_json, "atr", &atr, &m);
  status += extract_json_double_array(test_data_json, "sigx", &sigx, &m);
  status += extract_json_double_array(test_data_json, "w1p", &w1p, &m);
  status += extract_json_double_array(test_data_json, "dx", &dx_exp, &m);
  // status +=extract_json_double_array(test_data_json, "dx0", &dx0, &m);

  status += extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &n);

  status += extract_json_double(test_data_json, "fe", &fe);
  status += extract_json_double(test_data_json, "cgtol", &cgtol);
  status += extract_json_double(test_data_json, "tau", &tau);
  status += extract_json_int(test_data_json, "cgmaxiter", &cg_maxiter);

  l1c_dct1_setup(n, m, pix_idx, &ax_funs);

  h11p_data.one_by_fe = 1.0 / fe;
  h11p_data.one_by_fe_sqrd = 1.0 / (fe * fe);
  h11p_data.atr = atr;
  h11p_data.sigx = sigx;
  h11p_data.Dwork_1m = l1c_malloc_double(m);
  h11p_data.AtAx = ax_funs.AtAx;

  DWORK4 = l1c_malloc_double_2D(4, m);
  dx0 = l1c_malloc_double(m);
  l1c_init_vec(m, dx0, 0);
  double* dx_by_nrm = l1c_malloc_double(4 * m);
  double* dx_by_nrm_exp = l1c_malloc_double(4 * m);
  if (status || !DWORK4 || !dx_by_nrm || !dx_by_nrm_exp) {
    fprintf(stderr, "error allocating memory or reading JSON data in test_cgsolve_h11p\n");
    status += 1;
    goto exit;
  }

  cgp.max_iter = cg_maxiter;
  cgp.tol = cgtol;
  cgp.verbose = 0;

  l1c_cgsolve(m, dx0, w1p, DWORK4, _l1c_l1qc_H11pfun, &h11p_data, &cgr, cgp);

  double nrm_exp = cblas_dnrm2(m, dx_exp, 1);
  double nrm = cblas_dnrm2(m, dx0, 1);

  for (int i = 0; i < m; i++) {
    dx_by_nrm_exp[i] = dx_exp[i] / nrm_exp;
    dx_by_nrm[i] = dx0[i] / nrm;
  }
  ck_assert_double_array_eq_tol(m, dx_by_nrm_exp, dx_by_nrm, TOL_DOUBLE * 10);

exit:
  l1c_free_double(dx_by_nrm_exp);
  l1c_free_double(dx_by_nrm);
  l1c_free_double(atr);
  l1c_free_double(sigx);
  l1c_free_double(w1p);
  l1c_free_double(dx_exp);
  free(pix_idx);
  l1c_free_double(dx0);
  l1c_free_double_2D(4, DWORK4);
  l1c_free_double(h11p_data.Dwork_1m);
  ax_funs.destroy();

  free(fpath);

  cJSON_Delete(test_data_json);

  if (status) {
    ck_abort();
  }
}
END_TEST

/* Test the matrix multiplication function, Ax_sym.*/
START_TEST(test_cgsolve_Ax_sym) {

  char* fpath = fullfile(test_data_dir, "ax_sym.json");

  double *A, *x, *y_exp, *y;

  l1c_int N = 0, na = 0, status = 0;

  if (load_file_to_json(fpath, &test_data_json)) {
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  // Inputs to get_gradient
  status += extract_json_double_array(test_data_json, "A", &A, &na);
  status += extract_json_double_array(test_data_json, "x", &x, &N);
  status += extract_json_double_array(test_data_json, "y", &y_exp, &N);

  y = l1c_malloc_double(N);
  if ((!y) | status) {
    fprintf(stderr, "Error allocating memory\n");
    status += 1;
    goto exit;
  }

  Ax_sym(N, x, y, A);

  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE_SUPER * 10);

  // Now, y should be non-zero. Should get the same answer. Regression against having beta !=0,
  // because dspmv computes alpha*A*x + b*y
  Ax_sym(N, x, y, A);
  ck_assert_double_array_eq_tol(N, y_exp, y, TOL_DOUBLE_SUPER * 10);

exit:
  l1c_free_double(A);
  l1c_free_double(x);
  l1c_free_double(y_exp);
  l1c_free_double(y);

  free(fpath);
  cJSON_Delete(test_data_json);

  if (status) {
    ck_abort();
  }
}
END_TEST

/* ----------------------------- cgsolve_diag_precond--------------------------------*/
START_TEST(test_cgsolve_diag_precond) {
  double tol = 0.0; //= 1e-6;
  l1c_int max_iter;
  l1c_CgParams cgp;
  l1c_CgResults cgr;

  double *A = NULL, *x = NULL, *x_exp = NULL, *b = NULL;
  double **Dwork = NULL, *M_inv_diag = NULL;
  int status = 0;
  l1c_int N = 0, na = 0;
  if (load_small_data(&A, &x_exp, &b, &N, &na, &max_iter, &tol)) {
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = l1c_malloc_double(N);
  M_inv_diag = l1c_malloc_double(N);
  Dwork = l1c_malloc_double_2D(5, N);
  if (!Dwork || !x || !M_inv_diag) {
    status += 1;
    goto exit;
  }

  /* Must initialize x now, with warm starting. */
  for (int i = 0; i < N; i++) {
    x[i] = 0.0;
    M_inv_diag[i] = 1.0;
  }
  l1c_cgsolve_diag_precond(N, x, b, M_inv_diag, Dwork, Ax_sym, A, &cgr, cgp);

  ck_assert_double_array_eq_tol(N, x_exp, x, TOL_DOUBLE);
  goto exit;

exit:
  l1c_free_double(A);
  l1c_free_double(x);
  l1c_free_double(x_exp);
  l1c_free_double(b);
  l1c_free_double_2D(5, Dwork);
  l1c_free_double(M_inv_diag);

  cJSON_Delete(test_data_json);

  if (status) {
    ck_abort();
  }
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* cgsolve_suite(void) {
  Suite* s;

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
