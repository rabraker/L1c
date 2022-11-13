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

#define TOL_LARGE_DCT TOL_DOUBLE_SUPER

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

  /* Transform is n by mtot. mtot = mrow*mcol.*/
  l1c_int mrow;
  l1c_int mcol;
  l1c_int mtot;
  l1c_int n;

  int setup_status;
  char* fpath;

} DctData;

/* Global variable for all test cases.
 */

static DctData* dctd;
static l1c_AxFuns ax_funs;

/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(char* fname) {
  cJSON* test_data_json;
  int setup_status = 0;

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
  if (load_file_to_json(dctd->fpath, &test_data_json)) {
    fprintf(stderr, "Error loading data in test_dct\n");
    ck_abort();
  }

  setup_status += extract_json_double_array(test_data_json, "x_in", &dctd->x_in, &dctd->mtot);
  setup_status += extract_json_double_array(test_data_json, "y_in", &dctd->y_in, &dctd->n);
  setup_status += extract_json_double_array(test_data_json, "z_in", &dctd->z_in, &dctd->mtot);

  setup_status +=
      extract_json_double_array(test_data_json, "MtEt_EMx", &dctd->MtEty_EMx_exp, &dctd->mtot);

  setup_status += extract_json_double_array(test_data_json, "MtEty", &dctd->MtEty_exp, &dctd->mtot);

  setup_status += extract_json_double_array(test_data_json, "EMx", &dctd->EMx_exp, &dctd->n);

  setup_status += extract_json_double_array(test_data_json, "Mx", &dctd->Mx_exp, &dctd->mtot);
  setup_status += extract_json_double_array(test_data_json, "Mty", &dctd->Mty_exp, &dctd->mtot);
  setup_status += extract_json_double_array(test_data_json, "Ex", &dctd->Ex_exp, &dctd->n);
  setup_status += extract_json_double_array(test_data_json, "Ety", &dctd->Ety_exp, &dctd->mtot);

  setup_status += extract_json_int_array(test_data_json, "pix_idx", &dctd->pix_idx, &dctd->n);
  setup_status += extract_json_int(test_data_json, "mrow", &dctd->mrow);
  setup_status += extract_json_int(test_data_json, "mcol", &dctd->mcol);

  if (setup_status) {
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(dctd->mrow * dctd->mcol, dctd->mtot);

  dctd->MtEty_EMx_act = l1c_malloc_double(dctd->mtot);
  dctd->MtEty_act = l1c_malloc_double(dctd->mtot);
  dctd->EMx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mty_act = l1c_malloc_double(dctd->mtot);
  dctd->Ex_act = l1c_calloc_double(dctd->n);
  dctd->Ety_act = l1c_malloc_double(dctd->mtot);

  l1c_dct2_setup(dctd->n, dctd->mrow, dctd->mcol, dctd->pix_idx, &ax_funs);

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
  ax_funs.destroy();
  free(dctd->fpath);

  free(dctd);
}

static void setup_square(void) { setup("dct2_small.json"); }

static void setup_tall(void) { setup("dct2_small_tall.json"); }

static void setup_wide(void) { setup("dct2_small_wide.json"); }

static void setup_pure_square(void) { setup("dct2_small_pure_dct.json"); }

START_TEST(test_dct2_MtEt_EMx) {
  /* Provided and freed by setup_small() and teardown_small()*/
  ax_funs.AtAx(dctd->x_in, dctd->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(
      dctd->mtot, dctd->MtEty_EMx_exp, dctd->MtEty_EMx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_MtEty) {
  ax_funs.Aty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->MtEty_exp, dctd->MtEty_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_EMx) {
  ax_funs.Ax(dctd->x_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->EMx_exp, dctd->EMx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Mx) {
  ax_funs.Mx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Mx_exp, dctd->Mx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Wz_synth) {
  /* Should be the identity, for synthesis.*/
  ax_funs.Wz(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->x_in, dctd->Mx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Wtx_synth) {
  /* Should be the identity, for synthesis.*/
  ax_funs.Wtx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->x_in, dctd->Mx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Mty) {
  ax_funs.Mty(dctd->z_in, dctd->Mty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Mty_exp, dctd->Mty_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Ex) {
  ax_funs.Rx(dctd->x_in, dctd->Ex_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Ex_exp, dctd->Ex_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_dct2_Ety) {
  ax_funs.Rty(dctd->y_in, dctd->Ety_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Ety_exp, dctd->Ety_act, TOL_DOUBLE_SUPER);
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* dct2_suite(void) {
  Suite* s;

  TCase *tc_dct2_square, *tc_dct2_pure, *tc_dct2_tall, *tc_dct2_wide;
  // TCase  *tc_dct2_square;
  s = suite_create("dct2");

  tc_dct2_square = tcase_create("dct2_small");
  tcase_add_checked_fixture(tc_dct2_square, setup_square, teardown);
  tcase_add_test(tc_dct2_square, test_dct2_MtEt_EMx);
  tcase_add_test(tc_dct2_square, test_dct2_MtEty);
  tcase_add_test(tc_dct2_square, test_dct2_EMx);
  tcase_add_test(tc_dct2_square, test_dct2_Mx);
  tcase_add_test(tc_dct2_square, test_dct2_Mty);
  tcase_add_test(tc_dct2_square, test_dct2_Ex);
  tcase_add_test(tc_dct2_square, test_dct2_Ety);
  tcase_add_test(tc_dct2_square, test_dct2_Wz_synth);
  tcase_add_test(tc_dct2_square, test_dct2_Wtx_synth);

  suite_add_tcase(s, tc_dct2_square);

  /* The test check what happens when pix_idx is all ones, ie,
     the identitiy. That is, these are equivalent to just taking the dct/idct
     with not sampling.
  */

  tc_dct2_pure = tcase_create("dct2_pure");
  tcase_add_checked_fixture(tc_dct2_pure, setup_pure_square, teardown);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEt_EMx);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEty);
  tcase_add_test(tc_dct2_pure, test_dct2_EMx);
  tcase_add_test(tc_dct2_pure, test_dct2_Mx);
  tcase_add_test(tc_dct2_pure, test_dct2_Mty);
  tcase_add_test(tc_dct2_pure, test_dct2_Ex);
  tcase_add_test(tc_dct2_pure, test_dct2_Ety);
  tcase_add_test(tc_dct2_pure, test_dct2_Wz_synth);
  tcase_add_test(tc_dct2_pure, test_dct2_Wtx_synth);

  suite_add_tcase(s, tc_dct2_pure);

  /* DCT2 of tall, skinny matrix.
   */

  tc_dct2_tall = tcase_create("dct2_tall");
  tcase_add_checked_fixture(tc_dct2_tall, setup_tall, teardown);
  tcase_add_test(tc_dct2_tall, test_dct2_MtEt_EMx);
  tcase_add_test(tc_dct2_tall, test_dct2_MtEty);
  tcase_add_test(tc_dct2_tall, test_dct2_EMx);
  tcase_add_test(tc_dct2_tall, test_dct2_Mx);
  tcase_add_test(tc_dct2_tall, test_dct2_Mty);
  tcase_add_test(tc_dct2_tall, test_dct2_Ex);
  tcase_add_test(tc_dct2_tall, test_dct2_Ety);
  tcase_add_test(tc_dct2_tall, test_dct2_Wz_synth);
  tcase_add_test(tc_dct2_tall, test_dct2_Wtx_synth);

  suite_add_tcase(s, tc_dct2_tall);

  /* DCT2 of wide matrix.
   */
  tc_dct2_wide = tcase_create("dct2_wide");
  tcase_add_checked_fixture(tc_dct2_wide, setup_wide, teardown);
  tcase_add_test(tc_dct2_wide, test_dct2_MtEt_EMx);
  tcase_add_test(tc_dct2_wide, test_dct2_MtEty);
  tcase_add_test(tc_dct2_wide, test_dct2_EMx);
  tcase_add_test(tc_dct2_wide, test_dct2_Mx);
  tcase_add_test(tc_dct2_wide, test_dct2_Mty);
  tcase_add_test(tc_dct2_wide, test_dct2_Ex);
  tcase_add_test(tc_dct2_wide, test_dct2_Ety);
  tcase_add_test(tc_dct2_wide, test_dct2_Wz_synth);
  tcase_add_test(tc_dct2_wide, test_dct2_Wtx_synth);

  suite_add_tcase(s, tc_dct2_wide);

  return s;
}
