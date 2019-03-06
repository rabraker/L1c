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

ImgDiffData *dd;


static void setup(ImgDiffData *ddat){
  cJSON *test_data_json;
  int setup_status=0;
  int tmp=0, NM=0;
  if (load_file_to_json(ddat->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in TV_suite\n");
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "A", &ddat->A, &ddat->NM);
  setup_status +=extract_json_double_array(test_data_json, "DxA", &ddat->DxA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyA", &ddat->DyA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DxTA", &ddat->DxTA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyTA", &ddat->DyTA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DxTDxA", &ddat->DxTDxA_exp, &tmp);
  setup_status +=extract_json_double_array(test_data_json, "DyTDyA", &ddat->DyTDyA_exp, &tmp);

  setup_status +=extract_json_double(test_data_json, "alpha", &ddat->alpha);
  setup_status +=extract_json_int(test_data_json, "M", &ddat->M);
  setup_status +=extract_json_int(test_data_json, "N", &ddat->N);

  if (setup_status){
    fprintf(stderr, "Error loading json into test data from file: \n %s. Aborting.\n",
            ddat->fpath);
    ck_abort();
  }

  NM = ddat->NM;

  ddat->DxA_act = malloc_double(NM);
  ddat->DyA_act = malloc_double(NM);
  ddat->DxTA_act = malloc_double(NM);
  ddat->DyTA_act = malloc_double(NM);

  ddat->DxTDxA_act = malloc_double(NM);
  ddat->DyTDyA_act = malloc_double(NM);

  if (!(ddat->DxA_act) || !(ddat->DyA_act) || !(ddat->DxTA_act)
      || !(ddat->DyTA_act) || !(ddat->DxTDxA_act) || !(ddat->DyTDyA_act)){
    fprintf(stderr, "Error allocating memory. Aborting.\n");
    ck_abort();
  }
  cJSON_Delete(test_data_json);
}

static void teardown(ImgDiffData *ddat){
  free_double(ddat->DxA_act);
  free_double(ddat->DyA_act);
  free_double(ddat->DxTA_act);
  free_double(ddat->DyTA_act);

  free_double(ddat->DxTDxA_act);
  free_double(ddat->DyTDyA_act);
}

void setup_square(void){
  dd = malloc(sizeof(ImgDiffData));

  char *fpath_dd = fullfile(test_data_dir, "TV_data_square_small.json");

  sprintf(dd->fpath, "%s", fpath_dd);
  setup(dd);
}

void teardown_square(void){
  teardown(dd);
  free(dd);
}

void setup_tall(void){
  dd = malloc(sizeof(ImgDiffData));

  char *fpath_dd = fullfile(test_data_dir, "TV_data_tall_skinny.json");

  sprintf(dd->fpath, "%s", fpath_dd);
  setup(dd);
}


void teardown_tall(void){
  teardown(dd);
  free(dd);
}


void setup_short(void){
  dd = malloc(sizeof(ImgDiffData));

  char *fpath_dd = fullfile(test_data_dir, "TV_data_short_wide.json");

  sprintf(dd->fpath, "%s", fpath_dd);
  setup(dd);
}

void teardown_short(void){
  teardown(dd);
  free(dd);
}


START_TEST(test_l1c_Dx){

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i=0; i<N*M; i++){
    dd->DxA_act[i] =1;
  }

  l1c_Dx(N, M, dd->A, dd->DxA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DxA_exp,
                                dd->DxA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST

START_TEST(test_l1c_DxT){

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i=0; i<N*M; i++){
    dd->DxTA_act[i] =1;
  }

  l1c_DxT(N, M, dd->alpha, dd->A, dd->DxTA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DxTA_exp,
                                dd->DxTA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST

START_TEST(test_l1c_DxTDx){

  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i=0; i<N*M; i++){
    dd->DxTDxA_act[i] =1;
  }

  l1c_DxTDx(N, M, dd->alpha, dd->A, dd->DxTDxA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DxTDxA_exp,
                                dd->DxTDxA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


START_TEST(test_l1c_Dy){
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  for (int i=0; i<N*M; i++){
    dd->DyA_act[i] =0;
  }

  l1c_Dy(N, M, dd->A, dd->DyA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DyA_exp,
                                dd->DyA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST



START_TEST(test_l1c_DyT){
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  l1c_DyT(N, M, dd->alpha, dd->A, dd->DyTA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DyTA_exp,
                                dd->DyTA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


START_TEST(test_l1c_DyTDy){
  l1c_int N = dd->N;
  l1c_int M = dd->M;

  // double lap_y[12];
  l1c_DyTDy(N, M, dd->alpha, dd->A, dd->DyTDyA_act);

  ck_assert_double_array_eq_tol(N*M, dd->DyTDyA_exp,
                                dd->DyTDyA_act, TOL_DOUBLE_SUPER*10);
}
END_TEST


/* Add all the test cases to our suite
 */
Suite *TV_suite(void)
{
  Suite *s;

  TCase *tc_square, *tc_tall, *tc_short;
  s = suite_create("TV");

  tc_square = tcase_create("TV_square");
  tcase_add_checked_fixture(tc_square, setup_square, teardown_square);

  tcase_add_test(tc_square, test_l1c_Dx);
  tcase_add_test(tc_square, test_l1c_DxT);
  tcase_add_test(tc_square, test_l1c_DxTDx);

  tcase_add_test(tc_square, test_l1c_Dy);
  tcase_add_test(tc_square, test_l1c_DyT);
  tcase_add_test(tc_square, test_l1c_DyTDy);
  suite_add_tcase(s, tc_square);

  tc_tall = tcase_create("TV_tall");
  tcase_add_checked_fixture(tc_tall, setup_tall, teardown_tall);

  tcase_add_test(tc_tall, test_l1c_Dx);
  tcase_add_test(tc_tall, test_l1c_DxT);
  tcase_add_test(tc_tall, test_l1c_DxTDx);

  tcase_add_test(tc_tall, test_l1c_Dy);
  tcase_add_test(tc_tall, test_l1c_DyT);
  tcase_add_test(tc_tall, test_l1c_DyTDy);
  suite_add_tcase(s, tc_tall);

  tc_short = tcase_create("TV_short");
  tcase_add_checked_fixture(tc_short, setup_short, teardown_short);

  tcase_add_test(tc_short, test_l1c_Dx);
  tcase_add_test(tc_short, test_l1c_DxT);
  tcase_add_test(tc_short, test_l1c_DxTDx);

  tcase_add_test(tc_short, test_l1c_Dy);
  tcase_add_test(tc_short, test_l1c_DyT);
  tcase_add_test(tc_short, test_l1c_DyTDy);
  suite_add_tcase(s, tc_short);

  return s;

}
