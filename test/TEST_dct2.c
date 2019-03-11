/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */
#include "config.h"


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <fftw3.h>

#include "dct2.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"
#include "check_utils.h"
#include "l1c_common.h"
#include "l1c_memory.h"

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

/* Global variables for each test case: There are three test cases. The first
   uses randomly generated data with size Nx=50, the second is large, and from
   a simulation image, the third is is also random with Nx=50, but pix_idx corresonds to
   the Nx x Nx identity.
 */

DctData *dct2_small_data;
DctData *dct2_large_data;
DctData *dct2_pure_data;



/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
static void setup(DctData *dctd){
  cJSON *test_data_json;

  int setup_status=0;
  if (load_file_to_json(dctd->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_dct\n");
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "x_in", &dctd->x_in, &dctd->Nx);
  setup_status +=extract_json_double_array(test_data_json, "y_in", &dctd->y_in, &dctd->Ny);

  setup_status +=extract_json_double_array(test_data_json, "MtEt_EMx", &dctd->MtEty_EMx_exp, &dctd->Nx);

  setup_status +=extract_json_double_array(test_data_json, "MtEty", &dctd->MtEty_exp, &dctd->Nx);

  setup_status +=extract_json_double_array(test_data_json, "EMx", &dctd->EMx_exp, &dctd->Ny);

  setup_status +=extract_json_int_array(test_data_json, "pix_idx", &dctd->pix_idx, &dctd->Npix);
  setup_status +=extract_json_int(test_data_json, "N", &dctd->N);
  setup_status +=extract_json_int(test_data_json, "M", &dctd->M);

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(dctd->N * dctd->M, dctd->Nx);

  dctd->MtEty_EMx_act = malloc_double(dctd->Nx);
  dctd->MtEty_act = malloc_double(dctd->Nx);
  dctd->EMx_act = malloc_double(dctd->Nx);

  dct2_setup(dctd->N, dctd->M, dctd->Npix, dctd->pix_idx);

  cJSON_Delete(test_data_json);
}

static void teardown(DctData *dctd){
  free_double(dctd->x_in);
  free_double(dctd->y_in);

  free_double(dctd->MtEty_EMx_exp);
  free_double(dctd->MtEty_exp);
  free_double(dctd->EMx_exp);

  free_double(dctd->MtEty_EMx_act);
  free_double(dctd->MtEty_act);
  free_double(dctd->EMx_act);


  free(dctd->pix_idx);
  dct2_destroy();


#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}

static void setup_small(void){
  dct2_small_data = malloc(sizeof(DctData));

  char *fpath_dct2_small = fullfile(test_data_dir, "dct2_small.json");

  sprintf(dct2_small_data->fpath, "%s", fpath_dct2_small);
  setup(dct2_small_data);

  free(fpath_dct2_small);
}

static void teardown_small(void){
  teardown(dct2_small_data);
  free(dct2_small_data);
}



static void setup_pure_dct(void){
  dct2_pure_data = malloc(sizeof(DctData));

  char *fpath_dct2_pure = fullfile(test_data_dir, "dct2_small_pure_dct.json");
  sprintf(dct2_pure_data->fpath, "%s", fpath_dct2_pure);
  setup(dct2_pure_data);

  free(fpath_dct2_pure);
}

static void teardown_pure(void){
  teardown(dct2_pure_data);
  free(dct2_pure_data);

}


START_TEST(test_dct2_MtEt_EMx_small)
{
  /* Provided and freed by setup_small() and teardown_small()*/
  dct2_MtEt_EMx(dct2_small_data->x_in, dct2_small_data->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dct2_small_data->Nx, dct2_small_data->MtEty_EMx_exp,
                                dct2_small_data->MtEty_EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST(test_dct2_MtEty_small)
{
  dct2_MtEty(dct2_small_data->y_in, dct2_small_data->MtEty_act);

  ck_assert_double_array_eq_tol(dct2_small_data->Nx, dct2_small_data->MtEty_exp,
                                dct2_small_data->MtEty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_dct2_EMx_new_small)
{
  dct2_EMx(dct2_small_data->x_in, dct2_small_data->EMx_act);

  ck_assert_double_array_eq_tol(dct2_small_data->Ny, dct2_small_data->EMx_exp,
                                dct2_small_data->EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


/* The next three tests check what happens when pix_idx is all ones, ie,
 the identitiy. That is, these are equivalent to just taking the dct/idct
with not sampling.
*/
START_TEST(test_dct2_MtEt_EMx_pure)
{
   /* Test the multiplication (EM)^T * (E*M) * x */
  dct2_MtEt_EMx(dct2_pure_data->x_in, dct2_pure_data->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dct2_pure_data->Nx, dct2_pure_data->MtEty_EMx_exp,
                                dct2_pure_data->MtEty_EMx_act, TOL_DOUBLE_SUPER);


}
END_TEST


START_TEST(test_dct2_MtEty_pure)
{
  dct2_MtEty(dct2_pure_data->y_in, dct2_pure_data->MtEty_act);

  ck_assert_double_array_eq_tol(dct2_pure_data->Nx, dct2_pure_data->MtEty_exp,
                                dct2_pure_data->MtEty_act, TOL_DOUBLE_SUPER);


}
END_TEST




START_TEST(test_dct2_EMx_new_pure)
{
  dct2_EMx(dct2_pure_data->x_in, dct2_pure_data->EMx_act);

  ck_assert_double_array_eq_tol(dct2_pure_data->Ny, dct2_pure_data->EMx_exp,
                                dct2_pure_data->EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST




/* Add all the test cases to our suite
 */
Suite *dct2_suite(void)
{
  Suite *s;

  TCase *tc_dct2_small, *tc_dct2_pure;
  // TCase  *tc_dct2_small;
  s = suite_create("dct2");


  tc_dct2_small = tcase_create("dct2_small");
  tcase_add_checked_fixture(tc_dct2_small, setup_small, teardown_small);
  tcase_add_test(tc_dct2_small, test_dct2_MtEt_EMx_small);
  tcase_add_test(tc_dct2_small, test_dct2_MtEty_small);
  tcase_add_test(tc_dct2_small, test_dct2_EMx_new_small);

  suite_add_tcase(s, tc_dct2_small);


  tc_dct2_pure = tcase_create("dct2_pure");
  tcase_add_checked_fixture(tc_dct2_pure, setup_pure_dct, teardown_pure);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEt_EMx_pure);
  tcase_add_test(tc_dct2_pure, test_dct2_EMx_new_pure);
  tcase_add_test(tc_dct2_pure, test_dct2_MtEty_pure);

  suite_add_tcase(s, tc_dct2_pure);




  return s;

}
