/*
  Tests for the Bregman Splitting optimizations.

 */
#include "config.h"
#define CK_FLOATING_DIG 20

#include "cblas.h"
#include <check.h>
#include <cjson/cJSON.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Tolerances and things */
#include "bregman.h"
#include "check_utils.h"
#include "json_utils.h"
#include "l1c.h"
#include "l1c_math.h"
#include "test_constants.h"

extern char* fullfile(char* base_path, char* name);
extern char* test_data_dir;

typedef struct BregData {
  l1c_int NM;
  l1c_int n;
  l1c_int m;
  double* x;
  double* y;
  double* z;
  double* x_shrunk_exp;
  double* dx;
  double* dy;
  double* bx;
  double* by;
  double* f;
  double* rhs_exp;
  double* b;
  double* Hsolveb;
  double* Hessx;
  double* H_diag_exp;
  double gamma;
  double mu;
  double lam;
  char* fpath;
} BregData;

static BregData* BD;

static void setup(void) {
  BD = malloc(sizeof(BregData));
  BD->fpath = fullfile(test_data_dir, "bregman.json");
  if (!BD->fpath)
    ck_abort();
  cJSON* td_json;
  int setup_status = 0;
  int tmp = 0;
  if (load_file_to_json(BD->fpath, &td_json)) {
    fprintf(stderr, "Error loading data in %s (from %s)\n", __func__, __FILE__);
    ck_abort();
  }

  setup_status += extract_json_double_array(td_json, "x", &BD->x, &BD->NM);
  setup_status += extract_json_double_array(td_json, "y", &BD->y, &tmp);
  setup_status += extract_json_double_array(td_json, "z", &BD->z, &tmp);
  setup_status +=
      extract_json_double_array(td_json, "x_shrunk", &BD->x_shrunk_exp, &tmp);
  setup_status += extract_json_double_array(td_json, "dx", &BD->dx, &tmp);
  setup_status += extract_json_double_array(td_json, "dy", &BD->dy, &tmp);
  setup_status += extract_json_double_array(td_json, "bx", &BD->bx, &tmp);
  setup_status += extract_json_double_array(td_json, "by", &BD->by, &tmp);
  setup_status += extract_json_double_array(td_json, "f", &BD->f, &tmp);
  setup_status += extract_json_double_array(td_json, "rhs", &BD->rhs_exp, &tmp);
  setup_status += extract_json_double_array(td_json, "b", &BD->b, &tmp);
  setup_status += extract_json_double_array(td_json, "Hsolveb", &BD->Hsolveb, &tmp);
  setup_status += extract_json_double_array(td_json, "Hessx", &BD->Hessx, &tmp);
  setup_status += extract_json_double_array(td_json, "H_diag", &BD->H_diag_exp, &tmp);

  setup_status += extract_json_double(td_json, "gamma", &BD->gamma);
  setup_status += extract_json_double(td_json, "mu", &BD->mu);
  setup_status += extract_json_double(td_json, "lam", &BD->lam);
  setup_status += extract_json_int(td_json, "m", &BD->m);
  setup_status += extract_json_int(td_json, "n", &BD->n);

  if (setup_status) {
    fprintf(stderr,
            "Error loading json into test data from file: \n %s. In %s, Aborting.\n",
            BD->fpath,
            __func__);
    ck_abort();
  }

  cJSON_Delete(td_json);
}

static void teardown(void) {
  l1c_free_double(BD->x);
  l1c_free_double(BD->y);
  l1c_free_double(BD->z);
  l1c_free_double(BD->x_shrunk_exp);
  l1c_free_double(BD->dx);
  l1c_free_double(BD->dy);
  l1c_free_double(BD->bx);
  l1c_free_double(BD->by);
  l1c_free_double(BD->f);
  l1c_free_double(BD->rhs_exp);
  l1c_free_double(BD->b);
  l1c_free_double(BD->Hsolveb);
  l1c_free_double(BD->Hessx);
  l1c_free_double(BD->H_diag_exp);
  free(BD->fpath);
  free(BD);
}

START_TEST(test_breg_anis_jacobi) {
  BregFuncs bfuncs = breg_get_functions();

  double* uk = l1c_malloc_double(BD->NM);
  double* dwork = l1c_malloc_double(BD->NM);
  double* D = l1c_malloc_double(BD->NM);
  int i = 0;

  // Must initialize uk, becuase it used recursively in breg_anis_jacobi.
  for (i = 0; i < BD->NM; i++) {
    uk[i] = BD->b[i];
  }
  bfuncs.hess_inv_diag(BD->n, BD->m, BD->mu, BD->lam, D);

  for (i = 0; i < (BD->NM); i++) {
    bfuncs.breg_anis_jacobi(BD->n, BD->m, uk, dwork, BD->b, D, BD->lam);
  }

  ck_assert_double_array_eq_tol(BD->NM, BD->Hsolveb, uk, TOL_DOUBLE);
  l1c_free_double(uk);
  l1c_free_double(dwork);
  l1c_free_double(D);
}
END_TEST

START_TEST(test_breg_anis_guass_seidel) {
  BregFuncs bfuncs = breg_get_functions();

  double* uk = l1c_malloc_double(BD->NM);
  cblas_dcopy(BD->NM, BD->b, 1, uk, 1);

  for (int i = 0; i < (BD->NM); i++) {
    bfuncs.breg_anis_guass_seidel(BD->n, BD->m, uk, BD->b, BD->mu, BD->lam);
  }

  ck_assert_double_array_eq_tol(BD->NM, BD->Hsolveb, uk, TOL_DOUBLE);
  l1c_free_double(uk);
}
END_TEST

START_TEST(test_breg_shrink1) {

  BregFuncs bfuncs = breg_get_functions();

  double* x_shrunk = NULL;
  l1c_int N = 0;

  N = BD->NM;
  x_shrunk = l1c_malloc_double(N);
  if ((!x_shrunk)) {
    fprintf(stderr, "Error allocating memory in %s.\n", __func__);
    ck_abort();
  }
  for (int i = 0; i < N; i++) {
    x_shrunk[i] = 0;
  }

  bfuncs.breg_shrink1(N, BD->x, x_shrunk, BD->gamma);

  ck_assert_double_array_eq_tol(N, BD->x_shrunk_exp, x_shrunk, TOL_DOUBLE_SUPER * 10);

  l1c_free_double(x_shrunk);
}
END_TEST

START_TEST(test_breg_rhs) {
  BregFuncs bfuncs = breg_get_functions();
  double *dwork1, *dwork2, *rhs;

  rhs = l1c_malloc_double(BD->NM);
  dwork1 = l1c_malloc_double(BD->NM);
  dwork2 = l1c_malloc_double(BD->NM);
  if ((!rhs || !dwork1 || !dwork2)) {
    fprintf(stderr, "Error allocating memory in %s.\n", __func__);
    ck_abort();
  }
  bfuncs.breg_anis_rhs(BD->n,
                       BD->m,
                       BD->f,
                       BD->dx,
                       BD->bx,
                       BD->dy,
                       BD->by,
                       rhs,
                       BD->mu,
                       BD->lam,
                       dwork1,
                       dwork2);
  ck_assert_double_array_eq_tol(BD->NM, BD->rhs_exp, rhs, TOL_DOUBLE_SUPER * 10);

  l1c_free_double(rhs);
  l1c_free_double(dwork1);
  l1c_free_double(dwork2);
}
END_TEST

START_TEST(test_breg_hess_inv_diag) {
  BregFuncs bfuncs = breg_get_functions();
  double* D;

  D = l1c_malloc_double(BD->NM);
  if ((!D)) {
    fprintf(stderr, "Error allocating memory in %s.\n", __func__);
    ck_abort();
  }

  bfuncs.hess_inv_diag(BD->n, BD->m, BD->mu, BD->lam, D);

  ck_assert_double_array_eq_tol(BD->NM, BD->H_diag_exp, D, TOL_DOUBLE_SUPER * 10);

  l1c_free_double(D);
}
END_TEST

START_TEST(test_breg_mxpy_z) {

  BregFuncs bfuncs = breg_get_functions();
  double* z;
  l1c_int N = BD->NM;

  z = l1c_malloc_double(N);
  if (!z) {
    fprintf(stderr, "Error allocating memory in %s.\n", __func__);
    ck_abort();
  }
  for (int i = 0; i < N; i++) {
    z[i] = 1;
  }
  bfuncs.breg_mxpy_z(N, BD->x, BD->y, z);

  ck_assert_double_array_eq_tol(N, BD->z, z, TOL_DOUBLE_SUPER * 10);

  l1c_free_double(z);
}
END_TEST

/*
  This test doesnt do much (not sure what to check), but
  is here so we can run valgrind on the function.
 */
START_TEST(test_breg_anis_TV) {
  char* fpath = NULL;
  cJSON* img_json = NULL;

  int n = 0, m = 0, NM = 0;
  double *uk = NULL, *img_vec = NULL;
  double mu = 5, tol = 0.001;
  int max_iter = 100, max_jac_iter = 1, status = 0;

  fpath = fullfile(test_data_dir, "bregman_img.json");
  if (!fpath)
    ck_abort();

  if (load_file_to_json(fpath, &img_json)) {
    fprintf(stderr, "Error loading data in %s.\n", __func__);
    status = L1C_FILE_READ_FAILURE;
    goto exit;
  }

  status += extract_json_double_array(img_json, "img_vec", &img_vec, &NM);
  status += extract_json_int(img_json, "n", &n);
  status += extract_json_int(img_json, "m", &m);
  uk = l1c_malloc_double(n * m);

  if (status || !uk) {
    status = L1C_OUT_OF_MEMORY;
    fprintf(stderr, "Error allocating memory in %s.\n", __func__);
    goto exit;
  }

  l1c_init_vec(n * m, uk, 0);
  int ret = l1c_breg_anistropic_TV(n, m, uk, img_vec, mu, tol, max_iter, max_jac_iter);

  ck_assert_int_eq(ret, 0);

exit:
  l1c_free_double(uk);
  l1c_free_double(img_vec);
  free(fpath);
  cJSON_Delete(img_json);

  if (status)
    ck_abort();
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* bregman_suite(void) {
  Suite* s;

  TCase *tc_breg, *tc_breg_anisTV;
  s = suite_create("bregman");
  tc_breg = tcase_create("breg");
  tcase_add_checked_fixture(tc_breg, setup, teardown);

  tcase_add_test(tc_breg, test_breg_shrink1);
  tcase_add_test(tc_breg, test_breg_mxpy_z);
  tcase_add_test(tc_breg, test_breg_hess_inv_diag);
  tcase_add_test(tc_breg, test_breg_rhs);
  tcase_add_test(tc_breg, test_breg_anis_guass_seidel);
  tcase_add_test(tc_breg, test_breg_anis_jacobi);
  suite_add_tcase(s, tc_breg);

  tc_breg_anisTV = tcase_create("breg_anisTV");
  tcase_add_test(tc_breg_anisTV, test_breg_anis_TV);
  suite_add_tcase(s, tc_breg_anisTV);
  return s;
}
