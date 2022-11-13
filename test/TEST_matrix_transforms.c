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

#include "check_utils.h"
#include "json_utils.h"
#include "l1c.h"
#include "l1c_math.h"
#include "test_constants.h"
#include <cjson/cJSON.h>

#define TOL_LARGE_DCT TOL_DOUBLE_SUPER

extern char* fullfile(char* base_path, char* name);
extern char* test_data_dir;

typedef struct MXfmData {
  double* A;
  double* AtAx_exp;
  double* AtAx_act;
  double* Aty_act;
  double* Aty_exp;
  double* Ax_act;
  double* Ax_exp;
  double* x_in;
  double* y_in;

  l1c_int Nrow;
  l1c_int Mcol;
  l1c_int NM;

  int setup_status;
  char fpath[256];

} MXfmData;

/* Global variable for each all test cases.  */
static MXfmData* mxfmd;

static l1c_AxFuns ax_funs;

/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(MXfmData* mxfm_dat) {
  cJSON* test_data_json;

  int setup_status = 0;
  if (load_file_to_json(mxfm_dat->fpath, &test_data_json)) {
    fprintf(stderr,
            "Error loading data from file %s (in %s at %s)\n",
            mxfm_dat->fpath,
            __FILE__,
            __func__);
    ck_abort();
  }
  setup_status += extract_json_double_array(test_data_json, "A", &mxfm_dat->A, &mxfm_dat->NM);
  setup_status += extract_json_double_array(test_data_json, "x", &mxfm_dat->x_in, &mxfm_dat->Mcol);
  setup_status += extract_json_double_array(test_data_json, "y", &mxfm_dat->y_in, &mxfm_dat->Nrow);

  setup_status +=
      extract_json_double_array(test_data_json, "AtAx", &mxfm_dat->AtAx_exp, &mxfm_dat->Mcol);

  setup_status +=
      extract_json_double_array(test_data_json, "Aty", &mxfm_dat->Aty_exp, &mxfm_dat->Mcol);

  setup_status +=
      extract_json_double_array(test_data_json, "Ax", &mxfm_dat->Ax_exp, &mxfm_dat->Nrow);

  if (setup_status) {
    fprintf(stderr,
            "Error Loading json into program data in %s at %s'. Aborting\n",
            __FILE__,
            __func__);
    ck_abort();
  }
  mxfm_dat->AtAx_act = l1c_malloc_double(mxfm_dat->Mcol);
  mxfm_dat->Aty_act = l1c_malloc_double(mxfm_dat->Mcol);
  mxfm_dat->Ax_act = l1c_malloc_double(mxfm_dat->Nrow);

  /* cblas_dgemv works like:
     y = Ax + beta*y. We had a bug where beta=0, but
     initializing the result vectors to 0 covers that up.
   */
  l1c_init_vec(mxfmd->Mcol, mxfmd->AtAx_act, 1);
  l1c_init_vec(mxfmd->Mcol, mxfmd->Aty_act, 1);
  l1c_init_vec(mxfmd->Nrow, mxfmd->Ax_act, 1);

  l1c_setup_matrix_transforms(mxfm_dat->Nrow, mxfm_dat->Mcol, mxfm_dat->A, &ax_funs);

  cJSON_Delete(test_data_json);
}

static void teardown(MXfmData* mxfm_dat) {

  l1c_free_double(mxfm_dat->A);
  l1c_free_double(mxfm_dat->x_in);
  l1c_free_double(mxfm_dat->y_in);

  l1c_free_double(mxfm_dat->AtAx_exp);
  l1c_free_double(mxfm_dat->Aty_exp);
  l1c_free_double(mxfm_dat->Ax_exp);

  l1c_free_double(mxfm_dat->AtAx_act);
  l1c_free_double(mxfm_dat->Aty_act);
  l1c_free_double(mxfm_dat->Ax_act);

  ax_funs.destroy();
}

static void setup_small(void) {
  mxfmd = malloc(sizeof(MXfmData));

  char* fpath_dct_small = fullfile(test_data_dir, "matrix_xfm_small_data.json");
  if (!fpath_dct_small)
    ck_abort();

  sprintf(mxfmd->fpath, "%s", fpath_dct_small);
  setup(mxfmd);

  free(fpath_dct_small);
}

static void teardown_small(void) {
  teardown(mxfmd);
  free(mxfmd);
}

START_TEST(test_AtAx) {
  /* Provided and freed by setup_small() and teardown_small()*/
  ax_funs.AtAx(mxfmd->x_in, mxfmd->AtAx_act);

  ck_assert_double_array_eq_tol(mxfmd->Mcol, mxfmd->AtAx_exp, mxfmd->AtAx_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_Aty) {
  ax_funs.Aty(mxfmd->y_in, mxfmd->Aty_act);

  ck_assert_double_array_eq_tol(mxfmd->Mcol, mxfmd->Aty_exp, mxfmd->Aty_act, TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_Ax) {
  ax_funs.Ax(mxfmd->x_in, mxfmd->Ax_act);

  ck_assert_double_array_eq_tol(mxfmd->Nrow, mxfmd->Ax_exp, mxfmd->Ax_act, TOL_DOUBLE_SUPER);
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* matrix_transform_suite(void) {
  Suite* s;

  TCase* tc_matrix_xfm;
  // TCase  *tc_dct_small;
  s = suite_create("matrix_transforms");

  tc_matrix_xfm = tcase_create("matrix_xfm_small");
  tcase_add_checked_fixture(tc_matrix_xfm, setup_small, teardown_small);
  tcase_add_test(tc_matrix_xfm, test_AtAx);
  tcase_add_test(tc_matrix_xfm, test_Aty);
  tcase_add_test(tc_matrix_xfm, test_Ax);
  suite_add_tcase(s, tc_matrix_xfm);

  return s;
}
