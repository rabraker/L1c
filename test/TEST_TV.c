/*
  Tests for Image Gradients and Laplacian.
 */

#include "config.h"
#define CK_FLOATING_DIG 20

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>


/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"

#include "l1c_common.h"
#include "TV.h"
#include "check_utils.h"

extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;


typedef struct ImDiffData{
  double *A;
  double *DxA_exp;
  double *DyA_exp;

  double *DxTA_exp;
  double *DyTA_exp;

  double *DxTDxA_exp;
  double *DyTDyA_exp;

  double *DxA_act;
  double *DyA_act;

  double *DyTA_act;
  double *DxTA_act;

  double *DxTDxA_act;
  double *DyTDyA_act;

  double alpha;
  int N;
  int M;
  int NM;
  char fpath[256];
}ImgDiffData;

ImgDiffData *dd_square_small;



static void setup(ImgDiffData *dd){
  cJSON *test_data_json;
  int setup_status=0;
  int tmp=0, NM=0;
  if (load_file_to_json(dd->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in TV_suite\n");
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "A", &dd->A, &dd->NM);
  setup_status +=extract_json_double_array(test_data_json, "DxA", &dd->DxA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyA", &dd->DyA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DxTA", &dd->DxTA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyTA", &dd->DyTA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DxTDxA", &dd->DxTDxA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyTDyA", &dd->DyTDyA_exp, &tmp);

  setup_status +=extract_json_double(test_data_json, "alpha", &dd->alpha);
  setup_status +=extract_json_int(test_data_json, "M", &dd->N);
  setup_status +=extract_json_int(test_data_json, "N", &dd->M);

  if (setup_status){
    fprintf(stderr, "Error loading json into test data from file: \n %s. Aborting.\n",
            dd->fpath);
    ck_abort();
  }

  NM = dd->NM;

  dd->DxA_act = malloc_double(NM);
  dd->DyA_act = malloc_double(NM);
  dd->DxTA_act = malloc_double(NM);
  dd->DyTA_act = malloc_double(NM);

  dd->DxTDxA_act = malloc_double(NM);
  dd->DyTDyA_act = malloc_double(NM);

  if (!(dd->DxA_act) || !(dd->DyA_act) || !(dd->DxTA_act)
      || !(dd->DyTA_act) || !(dd->DxTDxA_act) || !(dd->DyTDyA_act)){
    fprintf(stderr, "Error allocating memory. Aborting.\n");
    ck_abort();
  }
  cJSON_Delete(test_data_json);
}

static void teardown(ImgDiffData *dd){
  free_double(dd->DxA_act);
  free_double(dd->DyA_act);
  free_double(dd->DxTA_act);
  free_double(dd->DyTA_act);

  free_double(dd->DxTDxA_act);
  free_double(dd->DyTDyA_act);
}

void setup_small(void){
  dd_square_small = malloc(sizeof(ImgDiffData));

  char *fpath_dd_square_small = fullfile(test_data_dir, "TV_data_square_small.json");

  sprintf(dd_square_small->fpath, "%s", fpath_dd_square_small);
  setup(dd_square_small);
}

void teardown_small(void){
  teardown(dd_square_small);
  free(dd_square_small);
}

START_TEST(test_l1c_Dx){

  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  for (int i=0; i<N*M; i++){
    dd_square_small->DxA_act[i] =0;
  }

  l1c_Dx(N, M, dd_square_small->A, dd_square_small->DxA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DxA_exp,
                                dd_square_small->DxA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST

START_TEST(test_l1c_DxT){

  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  for (int i=0; i<N*M; i++){
    dd_square_small->DxTA_act[i] =0;
  }

  l1c_DxT(N, M, dd_square_small->alpha, dd_square_small->A, dd_square_small->DxTA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DxTA_exp,
                                dd_square_small->DxTA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST

START_TEST(test_l1c_DxTDx){

  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  for (int i=0; i<N*M; i++){
    dd_square_small->DxTDxA_act[i] =0;
  }

  l1c_DxT(N, M, dd_square_small->alpha, dd_square_small->A, dd_square_small->DxTDxA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DxTDxA_exp,
                                dd_square_small->DxTDxA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


START_TEST(test_l1c_Dy){
  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  for (int i=0; i<N*M; i++){
    dd_square_small->DyA_act[i] =0;
  }

  l1c_Dy(N, M, dd_square_small->A, dd_square_small->DyA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DyA_exp,
                                dd_square_small->DyA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST



START_TEST(test_l1c_DyT){
  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  l1c_DyT(N, M, dd_square_small->alpha, dd_square_small->A, dd_square_small->DyTA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DyTA_exp,
                                dd_square_small->DyTA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


START_TEST(test_l1c_DyTDy){
  l1c_int N = dd_square_small->N;
  l1c_int M = dd_square_small->M;

  // double lap_y[12];
  l1c_DyTDy(N, M, dd_square_small->alpha, dd_square_small->A, dd_square_small->DyTDyA_act);

  ck_assert_double_array_eq_tol(N*M, dd_square_small->DyTDyA_exp,
                                dd_square_small->DyTDyA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


/* Add all the test cases to our suite
 */
Suite *TV_suite(void)
{
  Suite *s;

  TCase *tc_diff;
  s = suite_create("TV");

  tc_diff = tcase_create("diff");
  tcase_add_checked_fixture(tc_diff, setup_small, teardown_small);

  tcase_add_test(tc_diff, test_l1c_Dx);
  tcase_add_test(tc_diff, test_l1c_DxT);
  tcase_add_test(tc_diff, test_l1c_DxTDx);

  tcase_add_test(tc_diff, test_l1c_Dy);
  tcase_add_test(tc_diff, test_l1c_DyT);
  tcase_add_test(tc_diff, test_l1c_DyTDy);

  suite_add_tcase(s, tc_diff);

  return s;

}
