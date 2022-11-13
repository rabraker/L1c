/*
  Tests for Image Gradients and Laplacian.
 */

#include "config.h"
#define CK_FLOATING_DIG 20

#include <check.h>
#include <math.h> //Constants
#include <stdio.h>
#include <stdlib.h>

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "json_utils.h"
#include <cjson/cJSON.h>

#include "TV.h"
#include "check_utils.h"
#include "l1c.h"

extern char* fullfile(char* base_path, char* name);
extern char* test_data_dir;

typedef struct ImDiffData {
  double* A;
  double* DxA_exp;
  double* DyA_exp;

  double* DxTA_exp;
  double* DyTA_exp;

  double* DxTDxA_exp;
  double* DyTDyA_exp;

  double* DxA_act;
  double* DyA_act;

  double* DyTA_act;
  double* DxTA_act;

  double* DxTDxA_act;
  double* DyTDyA_act;

  double alpha;
  int N;
  int M;
  int NM;

} ImgDiffData;

static ImgDiffData* dd;

static void setup(char* fname) {
  cJSON* test_data_json;
  int setup_status = 0;
  int tmp = 0, NM = 0;
  dd = malloc(sizeof(ImgDiffData));
  if (!dd) {
    fprintf(stderr, "Memory allocation failed in %s() of file %s\n", __func__, __FILE__);
    ck_abort();
  }
  char* fpath = fullfile(test_data_dir, fname);
  /*TODO: properly free prior allocations if these fail.*/
  if (!fpath) {
    fprintf(stderr, "Memory allocation failed in %s() of file %s\n", __func__, __FILE__);
    ck_abort();
  }
  if (load_file_to_json(fpath, &test_data_json)) {
    fprintf(stderr, "Error loading data in TV_suite\n");
    ck_abort();
  }

  setup_status += extract_json_double_array(test_data_json, "A", &dd->A, &dd->NM);
  setup_status += extract_json_double_array(test_data_json, "DxA", &dd->DxA_exp, &tmp);
  setup_status += extract_json_double_array(test_data_json, "DyA", &dd->DyA_exp, &tmp);
  setup_status += extract_json_double_array(test_data_json, "DxTA", &dd->DxTA_exp, &tmp);
  setup_status += extract_json_double_array(test_data_json, "DyTA", &dd->DyTA_exp, &tmp);
  setup_status += extract_json_double_array(test_data_json, "DxTDxA", &dd->DxTDxA_exp, &tmp);
  setup_status += extract_json_double_array(test_data_json, "DyTDyA", &dd->DyTDyA_exp, &tmp);

  setup_status += extract_json_double(test_data_json, "alpha", &dd->alpha);
  setup_status += extract_json_int(test_data_json, "M", &dd->M);
  setup_status += extract_json_int(test_data_json, "N", &dd->N);

  if (setup_status) {
    fprintf(stderr, "Error loading json test data from file: \n %s. Aborting.\n", fpath);
    ck_abort();
  }

  NM = dd->NM;

  dd->DxA_act = l1c_malloc_double(NM);
  dd->DyA_act = l1c_malloc_double(NM);
  dd->DxTA_act = l1c_malloc_double(NM);
  dd->DyTA_act = l1c_malloc_double(NM);

  dd->DxTDxA_act = l1c_malloc_double(NM);
  dd->DyTDyA_act = l1c_malloc_double(NM);

  if (!(dd->DxA_act) || !(dd->DyA_act) || !(dd->DxTA_act) || !(dd->DyTA_act) || !(dd->DxTDxA_act) ||
      !(dd->DyTDyA_act)) {
    fprintf(stderr, "Error allocating memory. Aborting.\n");
    ck_abort();
  }
  cJSON_Delete(test_data_json);
  free(fpath);
}

static void teardown(void) {
  l1c_free_double(dd->A);
  l1c_free_double(dd->DxA_exp);
  l1c_free_double(dd->DyA_exp);
  l1c_free_double(dd->DxTA_exp);
  l1c_free_double(dd->DyTA_exp);
  l1c_free_double(dd->DxTDxA_exp);
  l1c_free_double(dd->DyTDyA_exp);

  l1c_free_double(dd->DxA_act);
  l1c_free_double(dd->DyA_act);
  l1c_free_double(dd->DxTA_act);
  l1c_free_double(dd->DyTA_act);

  l1c_free_double(dd->DxTDxA_act);
  l1c_free_double(dd->DyTDyA_act);

  free(dd);
}

static void setup_square(void) { setup("TV_data_square_small.json"); }

static void setup_tall(void) { setup("TV_data_tall_skinny.json"); }

static void setup_short(void) { setup("TV_data_short_wide.json"); }

START_TEST(test_l1c_Dx) {

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i = 0; i < N * M; i++) {
    dd->DxA_act[i] = 1;
  }

  l1c_Dx(N, M, dd->alpha, dd->A, dd->DxA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DxA_exp, dd->DxA_act, TOL_DOUBLE_SUPER * 10);
}
END_TEST

START_TEST(test_l1c_DxT) {

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i = 0; i < N * M; i++) {
    dd->DxTA_act[i] = 1;
  }

  l1c_DxT(N, M, dd->alpha, dd->A, dd->DxTA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DxTA_exp, dd->DxTA_act, TOL_DOUBLE_SUPER * 10);
}
END_TEST

START_TEST(test_l1c_DxTDx) {

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i = 0; i < N * M; i++) {
    dd->DxTDxA_act[i] = 1;
  }

  l1c_DxTDx(N, M, dd->alpha, dd->A, dd->DxTDxA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DxTDxA_exp, dd->DxTDxA_act, TOL_DOUBLE_SUPER * 100);
}
END_TEST

START_TEST(test_l1c_Dy) {
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i = 0; i < N * M; i++) {
    dd->DyA_act[i] = 0;
  }

  l1c_Dy(N, M, dd->alpha, dd->A, dd->DyA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DyA_exp, dd->DyA_act, TOL_DOUBLE_SUPER * 10);
}
END_TEST

START_TEST(test_l1c_DyT) {
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  l1c_DyT(N, M, dd->alpha, dd->A, dd->DyTA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DyTA_exp, dd->DyTA_act, TOL_DOUBLE_SUPER * 10);
}
END_TEST

START_TEST(test_l1c_DyTDy) {
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  // double lap_y[12];
  l1c_DyTDy(N, M, dd->alpha, dd->A, dd->DyTDyA_act);

  ck_assert_double_array_eq_tol(N * M, dd->DyTDyA_exp, dd->DyTDyA_act, TOL_DOUBLE_SUPER * 100);
}
END_TEST

/* Add all the test cases to our suite
 */
Suite* TV_suite(void) {
  Suite* s;

  TCase *tc_square, *tc_tall, *tc_short;
  s = suite_create("TV");

  tc_square = tcase_create("TV_square");
  tcase_add_checked_fixture(tc_square, setup_square, teardown);

  tcase_add_test(tc_square, test_l1c_Dx);
  tcase_add_test(tc_square, test_l1c_DxT);
  tcase_add_test(tc_square, test_l1c_DxTDx);

  tcase_add_test(tc_square, test_l1c_Dy);
  tcase_add_test(tc_square, test_l1c_DyT);
  tcase_add_test(tc_square, test_l1c_DyTDy);
  suite_add_tcase(s, tc_square);

  tc_tall = tcase_create("TV_tall");
  tcase_add_checked_fixture(tc_tall, setup_tall, teardown);

  tcase_add_test(tc_tall, test_l1c_Dx);
  tcase_add_test(tc_tall, test_l1c_DxT);
  tcase_add_test(tc_tall, test_l1c_DxTDx);

  tcase_add_test(tc_tall, test_l1c_Dy);
  tcase_add_test(tc_tall, test_l1c_DyT);
  tcase_add_test(tc_tall, test_l1c_DyTDy);
  suite_add_tcase(s, tc_tall);

  tc_short = tcase_create("TV_short");
  tcase_add_checked_fixture(tc_short, setup_short, teardown);

  tcase_add_test(tc_short, test_l1c_Dx);
  tcase_add_test(tc_short, test_l1c_DxT);
  tcase_add_test(tc_short, test_l1c_DxTDx);

  tcase_add_test(tc_short, test_l1c_Dy);
  tcase_add_test(tc_short, test_l1c_DyT);
  tcase_add_test(tc_short, test_l1c_DyTDy);
  suite_add_tcase(s, tc_short);

  return s;
}
