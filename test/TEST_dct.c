/*
  Tests for the 1D DCT-based transforms.

 */
#include "config.h"


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <fftw3.h>

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"
#include "check_utils.h"
#include "l1c.h"

#define TOL_DCT TOL_DOUBLE_SUPER

extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;

typedef struct DctData{
  l1c_int *pix_idx;
  double *MtEty_EMx_exp;
  double *MtEty_EMx_act;
  double *MtEty_act;
  double *MtEty_exp;
  double *EMx_act;
  double *EMx_exp;

  double *Mx_act;
  double *Mx_exp;
  double *Mty_act;
  double *Mty_exp;

  double *Ex_act;
  double *Ex_exp;
  double *Ety_act;
  double *Ety_exp;

  double *x_in;
  double *y_in;
  double *z_in;

  l1c_int m;
  l1c_int n;

  int setup_status;
  char fpath[256];

}DctData;

/* Global variable for each all test cases.  */
DctData *dctd;

l1c_AxFuns ax_funs;


/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(DctData *dct_dat){
  cJSON *test_data_json;

  int setup_status=0;
  if (load_file_to_json(dct_dat->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_dct from file %s\n", dct_dat->fpath);
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "x_in", &dct_dat->x_in, &dct_dat->m);
  setup_status +=extract_json_double_array(test_data_json, "y_in", &dct_dat->y_in, &dct_dat->n);
  setup_status +=extract_json_double_array(test_data_json, "z_in", &dct_dat->z_in, &dct_dat->m);

  setup_status +=extract_json_double_array(test_data_json, "MtEt_EMx", &dct_dat->MtEty_EMx_exp, &dct_dat->m);

  setup_status +=extract_json_double_array(test_data_json, "MtEty", &dct_dat->MtEty_exp, &dct_dat->m);

  setup_status +=extract_json_double_array(test_data_json, "EMx", &dct_dat->EMx_exp, &dct_dat->n);

  setup_status +=extract_json_double_array(test_data_json, "Mx", &dct_dat->Mx_exp, &dct_dat->m);
  setup_status +=extract_json_double_array(test_data_json, "Mty", &dct_dat->Mty_exp, &dct_dat->m);
  setup_status +=extract_json_double_array(test_data_json, "Ex", &dct_dat->Ex_exp, &dct_dat->n);
  setup_status +=extract_json_double_array(test_data_json, "Ety", &dct_dat->Ety_exp, &dct_dat->m);

  setup_status +=extract_json_int_array(test_data_json, "pix_idx", &dct_dat->pix_idx, &dct_dat->n);

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  dct_dat->MtEty_EMx_act = l1c_malloc_double(dct_dat->m);
  dct_dat->MtEty_act = l1c_malloc_double(dct_dat->m);
  dct_dat->EMx_act = l1c_malloc_double(dct_dat->m);
  dct_dat->Mx_act = l1c_malloc_double(dct_dat->m);
  dct_dat->Mty_act = l1c_malloc_double(dct_dat->m);
  dct_dat->Ex_act = l1c_calloc_double(dct_dat->n);
  dct_dat->Ety_act = l1c_malloc_double(dct_dat->m);

  l1c_dct1_setup(dct_dat->m, dct_dat->n, dct_dat->pix_idx, &ax_funs);

  cJSON_Delete(test_data_json);
}


static void teardown(DctData *dct_dat){
  l1c_free_double(dct_dat->x_in);
  l1c_free_double(dct_dat->y_in);
  l1c_free_double(dct_dat->z_in);

  l1c_free_double(dct_dat->MtEty_EMx_exp);
  l1c_free_double(dct_dat->MtEty_exp);
  l1c_free_double(dct_dat->EMx_exp);
  l1c_free_double(dct_dat->Mx_exp);
  l1c_free_double(dct_dat->Mty_exp);
  l1c_free_double(dct_dat->Ex_exp);
  l1c_free_double(dct_dat->Ety_exp);

  l1c_free_double(dct_dat->MtEty_EMx_act);
  l1c_free_double(dct_dat->MtEty_act);
  l1c_free_double(dct_dat->EMx_act);
  l1c_free_double(dct_dat->Mx_act);
  l1c_free_double(dct_dat->Mty_act);
  l1c_free_double(dct_dat->Ex_act);
  l1c_free_double(dct_dat->Ety_act);


  free(dct_dat->pix_idx);
  ax_funs.destroy();

}


static void setup_small(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct_small = fullfile(test_data_dir, "dct_small.json");

  sprintf(dctd->fpath, "%s", fpath_dct_small);
  setup(dctd);

  free(fpath_dct_small);
}


static void teardown_small(void){
  teardown(dctd);
  free(dctd);
}


/* */
static void setup_large(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct_large = fullfile(test_data_dir, "dct_large.json");
  sprintf(dctd->fpath, "%s", fpath_dct_large);
  setup(dctd);

  free(fpath_dct_large);

}


static void teardown_large(void){
  teardown(dctd);
  free(dctd);
}


static void setup_pure(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct_pure = fullfile(test_data_dir, "dct_small_pure_dct.json");
  sprintf(dctd->fpath, "%s", fpath_dct_pure);
  setup(dctd);

  free(fpath_dct_pure);
}


static void teardown_pure(void){
  teardown(dctd);
  free(dctd);

}


START_TEST(test_dct_MtEt_EMx)
{
  /* Provided and freed by setup_small() and teardown_small()*/
  ax_funs.AtAx(dctd->x_in, dctd->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->MtEty_EMx_exp,
                                dctd->MtEty_EMx_act, TOL_DCT);

}
END_TEST


START_TEST(test_dct_MtEty)
{
  ax_funs.Aty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->MtEty_exp,
                                dctd->MtEty_act, TOL_DCT);

}
END_TEST


START_TEST(test_dct_EMx)
{
  ax_funs.Ax(dctd->x_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->EMx_exp,
                                dctd->EMx_act, TOL_DCT);

}
END_TEST


START_TEST(test_dct_Mx)
{
  ax_funs.Mx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Mx_exp,
                                dctd->Mx_act, TOL_DCT);

}
END_TEST


START_TEST(test_dct_Mty)
{
  ax_funs.Mty(dctd->z_in, dctd->Mty_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Mty_exp,
                                dctd->Mty_act, 2*TOL_DCT);

}
END_TEST


START_TEST(test_dct_Ex)
{
  ax_funs.Ex(dctd->x_in, dctd->Ex_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Ex_exp,
                                dctd->Ex_act, TOL_DCT);

}
END_TEST

START_TEST(test_dct_Ety)
{
  ax_funs.Ety(dctd->y_in, dctd->Ety_act);

  ck_assert_double_array_eq_tol(dctd->m, dctd->Ety_exp,
                                dctd->Ety_act, TOL_DCT);

}
END_TEST

/* Add all the test cases to our suite
 */
Suite *dct_suite(void)
{
  Suite *s;

  TCase *tc_dct_small, *tc_dct_large, *tc_dct_pure;
  // TCase  *tc_dct_small;
  s = suite_create("dct");

  tc_dct_small = tcase_create("dct_small");
  tcase_add_checked_fixture(tc_dct_small, setup_small, teardown_small);
  tcase_add_test(tc_dct_small, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_small, test_dct_MtEty);
  tcase_add_test(tc_dct_small, test_dct_EMx);
  tcase_add_test(tc_dct_small, test_dct_Mx);
  tcase_add_test(tc_dct_small, test_dct_Mty);
  tcase_add_test(tc_dct_small, test_dct_Ex);
  tcase_add_test(tc_dct_small, test_dct_Ety);
  suite_add_tcase(s, tc_dct_small);

  tc_dct_large = tcase_create("dct_large");
  tcase_add_checked_fixture(tc_dct_large, setup_large, teardown_large);
  tcase_add_test(tc_dct_large,test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_large, test_dct_MtEty);
  tcase_add_test(tc_dct_large, test_dct_EMx);
  tcase_add_test(tc_dct_large, test_dct_Mx);
  tcase_add_test(tc_dct_large, test_dct_Mty);
  tcase_add_test(tc_dct_large, test_dct_Ex);
  tcase_add_test(tc_dct_large, test_dct_Ety);
  suite_add_tcase(s, tc_dct_large);

  tc_dct_pure = tcase_create("dct_pure");
  tcase_add_checked_fixture(tc_dct_pure, setup_pure, teardown_pure);
  tcase_add_test(tc_dct_pure, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_pure, test_dct_MtEty);
  tcase_add_test(tc_dct_pure, test_dct_EMx);
  tcase_add_test(tc_dct_pure, test_dct_Mx);
  tcase_add_test(tc_dct_pure, test_dct_Mty);
  tcase_add_test(tc_dct_pure, test_dct_Ex);
  tcase_add_test(tc_dct_pure, test_dct_Ety);
  suite_add_tcase(s, tc_dct_pure);

  return s;

}
