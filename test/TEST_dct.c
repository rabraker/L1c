/*
  Tests for the 1D DCT-based transforms.

 */
#include "config.h"

#define CK_FLOATING_DIG 20
#include <check.h>
#include <fftw3.h>
#include <math.h> //Constants
#include <stdio.h>
#include <stdlib.h>

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "check_utils.h"
#include "json_utils.h"
#include "l1c.h"
#include <cjson/cJSON.h>

#define TOL_DCT TOL_DOUBLE_SUPER

extern char* fullfile(char* base_path, char* name);
extern char* test_data_dir;

typedef struct DctData {
  l1c_int* pix_idx;
  double* MtEty_EMx_exp;
  double* MtEty_EMx_act;
  double* MtEty_act;
  double* MtEty_exp;
  double* EMx_act;
  double* EMx_exp;

  double* Mx_act;
  double* Mx_exp;
  double* Mty_act;
  double* Mty_exp;

  double* Ex_act;
  double* Ex_exp;
  double* Ety_act;
  double* Ety_exp;

  double* x_in;
  double* y_in;
  double* z_in;

  l1c_int m;
  l1c_int n;

  int setup_status;
  char* fpath;

} DctData;

/* Global variable for each all test cases.  */
static DctData* dctd;

static l1c_AxFuns ax_funs;

/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(char* fname) {
  cJSON* test_data_json;
  dctd = malloc(sizeof(DctData));
  if (!dctd) {
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }
  dctd->fpath = fullfile(test_data_dir, fname);

  if (!dctd->fpath) {
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }

  int setup_status = 0;
  if (load_file_to_json(dctd->fpath, &test_data_json)) {
    fprintf(stderr, "Error loading data in test_dct from file %s\n", dctd->fpath);
    ck_abort();
  }

  setup_status += extract_json_double_array(test_data_json, "x_in", &dctd->x_in, &dctd->m);
  setup_status += extract_json_double_array(test_data_json, "y_in", &dctd->y_in, &dctd->n);
  setup_status += extract_json_double_array(test_data_json, "z_in", &dctd->z_in, &dctd->m);

  setup_status +=
      extract_json_double_array(test_data_json, "MtEt_EMx", &dctd->MtEty_EMx_exp, &dctd->m);

  setup_status += extract_json_double_array(test_data_json, "MtEty", &dctd->MtEty_exp, &dctd->m);

  setup_status += extract_json_double_array(test_data_json, "EMx", &dctd->EMx_exp, &dctd->n);

  setup_status += extract_json_double_array(test_data_json, "Mx", &dctd->Mx_exp, &dctd->m);
  setup_status += extract_json_double_array(test_data_json, "Mty", &dctd->Mty_exp, &dctd->m);
  setup_status += extract_json_double_array(test_data_json, "Ex", &dctd->Ex_exp, &dctd->n);
  setup_status += extract_json_double_array(test_data_json, "Ety", &dctd->Ety_exp, &dctd->m);

  setup_status += extract_json_int_array(test_data_json, "pix_idx", &dctd->pix_idx, &dctd->n);

  if (setup_status) {
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  dctd->MtEty_EMx_act = l1c_malloc_double(dctd->m);
  dctd->MtEty_act = l1c_malloc_double(dctd->m);
  dctd->EMx_act = l1c_malloc_double(dctd->m);
  dctd->Mx_act = l1c_malloc_double(dctd->m);
  dctd->Mty_act = l1c_malloc_double(dctd->m);
  dctd->Ex_act = l1c_calloc_double(dctd->n);
  dctd->Ety_act = l1c_malloc_double(dctd->m);

  l1c_dct1_setup(dctd->n, dctd->m, dctd->pix_idx, &ax_funs);

  cJSON_Delete(test_data_json);
}

static void teardown(void) {
  l1c_free_double(dctd->x_in);
  l1c_free_double(dctd->y_in);
  l1c_free_double(dctd->z_in);

  l1c_free_double(dctd->MtEty_EMx_exp);
  l1c_free_double(dctd->MtEty_exp);
  l1c_free_double(dctd->EMx_exp);
  l1c_free_double(dctd->Mx_exp);
  l1c_free_double(dctd->Mty_exp);
  l1c_free_double(dctd->Ex_exp);
  l1c_free_double(dctd->Ety_exp);

  l1c_free_double(dctd->MtEty_EMx_act);
  l1c_free_double(dctd->MtEty_act);
  l1c_free_double(dctd->EMx_act);
  l1c_free_double(dctd->Mx_act);
  l1c_free_double(dctd->Mty_act);
  l1c_free_double(dctd->Ex_act);
  l1c_free_double(dctd->Ety_act);

  free(dctd->pix_idx);

  free(dctd->fpath);
  ax_funs.destroy();
  free(dctd);
}

static void setup_small(void) { setup("dct_small.json"); }

static void setup_large(void) { setup("dct_large.json"); }

static void setup_pure(void) { setup("dct_small_pure_dct.json"); }

START_TEST(test_dct_MtEt_EMx) {
  ax_funs.AtAx(dctd->x_in, dctd->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->MtEty_EMx_exp, dctd->MtEty_EMx_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_MtEty) {
  ax_funs.Aty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->MtEty_exp, dctd->MtEty_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_EMx) {
  ax_funs.Ax(dctd->x_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->EMx_exp, dctd->EMx_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_Mx) {
  ax_funs.Mx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Mx_exp, dctd->Mx_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_Wtx_synth) {
  /* Should be the identity.*/
  ax_funs.Wtx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->x_in, dctd->Mx_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_Wz_synth) {
  /* Should be the identity.*/
  ax_funs.Wz(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->x_in, dctd->Mx_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_Mty) {
  ax_funs.Mty(dctd->z_in, dctd->Mty_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Mty_exp, dctd->Mty_act, 2 * TOL_DCT);
}
END_TEST

START_TEST(test_dct_Rx) {
  ax_funs.Rx(dctd->x_in, dctd->Ex_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Ex_exp, dctd->Ex_act, TOL_DCT);
}
END_TEST

START_TEST(test_dct_Rty) {
  ax_funs.Rty(dctd->y_in, dctd->Ety_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Ety_exp, dctd->Ety_act, TOL_DCT);
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* dct_suite(void) {
  Suite* s;

  TCase *tc_dct_small, *tc_dct_large, *tc_dct_pure;
  // TCase  *tc_dct_small;
  s = suite_create("dct");

  tc_dct_small = tcase_create("dct_small");
  tcase_add_checked_fixture(tc_dct_small, setup_small, teardown);
  tcase_add_test(tc_dct_small, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_small, test_dct_MtEty);
  tcase_add_test(tc_dct_small, test_dct_EMx);
  tcase_add_test(tc_dct_small, test_dct_Mx);
  tcase_add_test(tc_dct_small, test_dct_Mty);
  tcase_add_test(tc_dct_small, test_dct_Rx);
  tcase_add_test(tc_dct_small, test_dct_Rty);
  tcase_add_test(tc_dct_small, test_dct_Wz_synth);
  tcase_add_test(tc_dct_small, test_dct_Wtx_synth);
  suite_add_tcase(s, tc_dct_small);

  tc_dct_large = tcase_create("dct_large");
  tcase_add_checked_fixture(tc_dct_large, setup_large, teardown);
  tcase_add_test(tc_dct_large, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_large, test_dct_MtEty);
  tcase_add_test(tc_dct_large, test_dct_EMx);
  tcase_add_test(tc_dct_large, test_dct_Mx);
  tcase_add_test(tc_dct_large, test_dct_Mty);
  tcase_add_test(tc_dct_large, test_dct_Rx);
  tcase_add_test(tc_dct_large, test_dct_Rty);
  tcase_add_test(tc_dct_large, test_dct_Wz_synth);
  tcase_add_test(tc_dct_large, test_dct_Wtx_synth);
  suite_add_tcase(s, tc_dct_large);

  tc_dct_pure = tcase_create("dct_pure");
  tcase_add_checked_fixture(tc_dct_pure, setup_pure, teardown);
  tcase_add_test(tc_dct_pure, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_pure, test_dct_MtEty);
  tcase_add_test(tc_dct_pure, test_dct_EMx);
  tcase_add_test(tc_dct_pure, test_dct_Mx);
  tcase_add_test(tc_dct_pure, test_dct_Mty);
  tcase_add_test(tc_dct_pure, test_dct_Rx);
  tcase_add_test(tc_dct_pure, test_dct_Rty);
  tcase_add_test(tc_dct_pure, test_dct_Wz_synth);
  tcase_add_test(tc_dct_pure, test_dct_Wtx_synth);
  suite_add_tcase(s, tc_dct_pure);

  return s;
}
