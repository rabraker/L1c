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
#include "l1c_common.h"
#include "l1c_memory.h"
#include "l1c_transforms.h"

#ifdef _USEMKL_
#define TOL_LARGE_DCT 1e-10
#else
#define TOL_LARGE_DCT TOL_DOUBLE_SUPER
#endif

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
  double *x_in;
  double *y_in;

  l1c_int N;
  l1c_int M;
  l1c_int Nx;
  l1c_int Npix;
  l1c_int Ny;

  int setup_status;
  char fpath[256];

}DctData;

/* Global variable for all test cases.
 */

DctData *dctd;
L1cAxFuns ax_funs;

/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(DctData *dct_dat){
  cJSON *test_data_json;

  int setup_status=0;
  if (load_file_to_json(dct_dat->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_dct\n");
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "x_in", &dct_dat->x_in, &dct_dat->Nx);
  setup_status +=extract_json_double_array(test_data_json, "y_in", &dct_dat->y_in, &dct_dat->Ny);

  setup_status +=extract_json_double_array(test_data_json, "MtEt_EMx", &dct_dat->MtEty_EMx_exp, &dct_dat->Nx);

  setup_status +=extract_json_double_array(test_data_json, "MtEty", &dct_dat->MtEty_exp, &dct_dat->Nx);

  setup_status +=extract_json_double_array(test_data_json, "EMx", &dct_dat->EMx_exp, &dct_dat->Ny);

  setup_status +=extract_json_int_array(test_data_json, "pix_idx", &dct_dat->pix_idx, &dct_dat->Npix);
  setup_status +=extract_json_int(test_data_json, "N", &dct_dat->N);
  setup_status +=extract_json_int(test_data_json, "M", &dct_dat->M);

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(dct_dat->N * dct_dat->M, dct_dat->Nx);

  dct_dat->MtEty_EMx_act = malloc_double(dct_dat->Nx);
  dct_dat->MtEty_act = malloc_double(dct_dat->Nx);
  dct_dat->EMx_act = malloc_double(dct_dat->Nx);

  dct2_setup(dct_dat->N, dct_dat->M, dct_dat->Npix, dct_dat->pix_idx, &ax_funs);

  cJSON_Delete(test_data_json);
}

static void teardown(DctData *dct_dat){
  free_double(dct_dat->x_in);
  free_double(dct_dat->y_in);

  free_double(dct_dat->MtEty_EMx_exp);
  free_double(dct_dat->MtEty_exp);
  free_double(dct_dat->EMx_exp);

  free_double(dct_dat->MtEty_EMx_act);
  free_double(dct_dat->MtEty_act);
  free_double(dct_dat->EMx_act);


  free(dct_dat->pix_idx);
  ax_funs.destroy();


#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}

static void setup_square(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct2_small = fullfile(test_data_dir, "dct2_small.json");

  sprintf(dctd->fpath, "%s", fpath_dct2_small);
  setup(dctd);

  free(fpath_dct2_small);
}

static void teardown_square(void){
  teardown(dctd);
  free(dctd);
}

static void setup_tall(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct2_small = fullfile(test_data_dir, "dct2_small_tall.json");

  sprintf(dctd->fpath, "%s", fpath_dct2_small);
  setup(dctd);

  free(fpath_dct2_small);
}

static void teardown_tall(void){
  teardown(dctd);
  free(dctd);
}


static void setup_wide(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct2_small = fullfile(test_data_dir, "dct2_small_wide.json");

  sprintf(dctd->fpath, "%s", fpath_dct2_small);
  setup(dctd);

  free(fpath_dct2_small);
}

static void teardown_wide(void){
  teardown(dctd);
  free(dctd);
}


static void setup_pure_square(void){
  dctd = malloc(sizeof(DctData));

  char *fpath_dct2_pure = fullfile(test_data_dir, "dct2_small_pure_dct.json");
  sprintf(dctd->fpath, "%s", fpath_dct2_pure);
  setup(dctd);

  free(fpath_dct2_pure);
}

static void teardown_pure_square(void){
  teardown(dctd);
  free(dctd);

}


START_TEST(test_dct2_MtEt_EMx_small)
{
  /* Provided and freed by setup_small() and teardown_small()*/
  ax_funs.AtAx(dctd->x_in, dctd->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dctd->Nx, dctd->MtEty_EMx_exp,
                                dctd->MtEty_EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST(test_dct2_MtEty_small)
{
  ax_funs.Aty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->Nx, dctd->MtEty_exp,
                                dctd->MtEty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_dct2_EMx_new_small)
{
  ax_funs.Ax(dctd->x_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->Ny, dctd->EMx_exp,
                                dctd->EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST



/* Add all the test cases to our suite
 */
Suite *dct2_suite(void)
{
  Suite *s;

  TCase *tc_dct2_square, *tc_dct2_pure, *tc_dct2_tall, *tc_dct2_wide;
  // TCase  *tc_dct2_square;
  s = suite_create("dct2");


  tc_dct2_square = tcase_create("dct2_small");
  tcase_add_checked_fixture(tc_dct2_square, setup_square, teardown_square);
  tcase_add_test(tc_dct2_square, test_dct2_MtEt_EMx_small);
  tcase_add_test(tc_dct2_square, test_dct2_MtEty_small);
  tcase_add_test(tc_dct2_square, test_dct2_EMx_new_small);

  suite_add_tcase(s, tc_dct2_square);

  /* The test check what happens when pix_idx is all ones, ie,
     the identitiy. That is, these are equivalent to just taking the dct/idct
     with not sampling.
  */

  tc_dct2_pure = tcase_create("dct2_pure");
  tcase_add_checked_fixture(tc_dct2_pure, setup_pure_square, teardown_pure_square);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEt_EMx_small);
  tcase_add_test(tc_dct2_pure, test_dct2_EMx_new_small);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEty_small);

  suite_add_tcase(s, tc_dct2_pure);

  /* DCT2 of tall, skinny matrix.
  */

  tc_dct2_tall = tcase_create("dct2_tall");
  tcase_add_checked_fixture(tc_dct2_tall, setup_tall, teardown_tall);
  tcase_add_test(tc_dct2_tall, test_dct2_MtEt_EMx_small);
  tcase_add_test(tc_dct2_tall, test_dct2_EMx_new_small);
  tcase_add_test(tc_dct2_tall, test_dct2_MtEty_small);

  suite_add_tcase(s, tc_dct2_tall);

  /* DCT2 of wide matrix.
  */
  tc_dct2_wide = tcase_create("dct2_wide");
  tcase_add_checked_fixture(tc_dct2_wide, setup_wide, teardown_wide);
  tcase_add_test(tc_dct2_wide, test_dct2_MtEt_EMx_small);
  tcase_add_test(tc_dct2_wide, test_dct2_EMx_new_small);
  tcase_add_test(tc_dct2_wide, test_dct2_MtEty_small);

  suite_add_tcase(s, tc_dct2_wide);

  return s;

}
