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

// #include "test_data.h"
#include "dct.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"
#include "check_utils.h"
#include "l1c_common.h"

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

DctData *dct_small_data;
DctData *dct_large_data;
DctData *dct_pure_data;



/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
void setup(DctData *dctd){
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

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  dctd->MtEty_EMx_act = malloc_double(dctd->Nx);
  dctd->MtEty_act = malloc_double(dctd->Nx);
  dctd->EMx_act = malloc_double(dctd->Nx);

  dct_setup(dctd->Nx, dctd->Npix, dctd->pix_idx);

  cJSON_Delete(test_data_json);
}

void teardown(DctData *dctd){
  free_double(dctd->x_in);
  free_double(dctd->y_in);

  free_double(dctd->MtEty_EMx_exp);
  free_double(dctd->MtEty_exp);
  free_double(dctd->EMx_exp);

  free_double(dctd->MtEty_EMx_act);
  free_double(dctd->MtEty_act);
  free_double(dctd->EMx_act);


  free(dctd->pix_idx);
  dct_destroy();


#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}

void setup_small(void){
  dct_small_data = malloc(sizeof(DctData));

  char *fpath_dct_small = fullfile(test_data_dir, "dct_small.json");

  sprintf(dct_small_data->fpath, "%s", fpath_dct_small);
  setup(dct_small_data);

  free(fpath_dct_small);
}

void teardown_small(void){
  teardown(dct_small_data);
  free(dct_small_data);
}

/* */
void setup_large(void){
  dct_large_data = malloc(sizeof(DctData));

  char *fpath_dct_large = fullfile(test_data_dir, "dct_large.json");
  sprintf(dct_large_data->fpath, "%s", fpath_dct_large);
  setup(dct_large_data);

  free(fpath_dct_large);

}

void teardown_large(void){
  teardown(dct_large_data);
  free(dct_large_data);
}

void setup_pure_dct(void){
  dct_pure_data = malloc(sizeof(DctData));

  char *fpath_dct_pure = fullfile(test_data_dir, "dct_small_pure_dct.json");
  sprintf(dct_pure_data->fpath, "%s", fpath_dct_pure);
  setup(dct_pure_data);

  free(fpath_dct_pure);
}

void teardown_pure(void){
  teardown(dct_pure_data);
  free(dct_pure_data);

}


START_TEST(test_dct_MtEt_EMx_small)
{
  /* Provided and freed by setup_small() and teardown_small()*/
  dct_MtEt_EMx_new(dct_small_data->x_in, dct_small_data->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dct_small_data->Nx, dct_small_data->MtEty_EMx_exp,
                                dct_small_data->MtEty_EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST(test_dct_MtEty_small)
{
  dct_MtEty(dct_small_data->y_in, dct_small_data->MtEty_act);

  ck_assert_double_array_eq_tol(dct_small_data->Nx, dct_small_data->MtEty_exp,
                                dct_small_data->MtEty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_dct_EMx_new_small)
{
  dct_EMx_new(dct_small_data->x_in, dct_small_data->EMx_act);

  ck_assert_double_array_eq_tol(dct_small_data->Ny, dct_small_data->EMx_exp,
                                dct_small_data->EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


/* -----------------------------------------------------------
 --------------- The large set, from an actual image ---------
 */

START_TEST(test_dct_MtEt_EMx_large)
{

  dct_MtEt_EMx_new(dct_large_data->x_in, dct_large_data->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dct_large_data->Nx, dct_large_data->MtEty_EMx_exp,
                                dct_large_data->MtEty_EMx_act,  TOL_LARGE_DCT);
}
END_TEST


START_TEST(test_dct_MtEty_large)
{
  dct_MtEty(dct_large_data->y_in, dct_large_data->MtEty_act);

  ck_assert_double_array_eq_tol(dct_large_data->Nx, dct_large_data->MtEty_exp,
                                dct_large_data->MtEty_act, TOL_LARGE_DCT);

}
END_TEST

START_TEST(test_dct_EMx_large)
{

  dct_EMx_new(dct_large_data->x_in, dct_large_data->EMx_act);

  ck_assert_double_array_eq_tol(dct_large_data->Ny, dct_large_data->EMx_exp,
                                dct_large_data->EMx_act, TOL_LARGE_DCT);

}
END_TEST



/* The next three tests check what happens when pix_idx is all ones, ie,
 the identitiy. That is, these are equivalent to just taking the dct/idct
with not sampling.
*/
START_TEST(test_dct_MtEt_EMx_pure)
{
   /* Test the multiplication (EM)^T * (E*M) * x */
  dct_MtEt_EMx_new(dct_pure_data->x_in, dct_pure_data->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dct_pure_data->Nx, dct_pure_data->MtEty_EMx_exp,
                                dct_pure_data->MtEty_EMx_act, TOL_DOUBLE_SUPER);


}
END_TEST


START_TEST(test_dct_MtEty_pure)
{
  dct_MtEty(dct_pure_data->y_in, dct_pure_data->MtEty_act);

  ck_assert_double_array_eq_tol(dct_pure_data->Nx, dct_pure_data->MtEty_exp,
                                dct_pure_data->MtEty_act, TOL_DOUBLE_SUPER);


}
END_TEST




START_TEST(test_dct_EMx_new_pure)
{
  dct_EMx_new(dct_pure_data->x_in, dct_pure_data->EMx_act);

  ck_assert_double_array_eq_tol(dct_pure_data->Ny, dct_pure_data->EMx_exp,
                                dct_pure_data->EMx_act, TOL_DOUBLE_SUPER);

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
  tcase_add_test(tc_dct_small, test_dct_MtEt_EMx_small);
  tcase_add_test(tc_dct_small, test_dct_MtEty_small);
  tcase_add_test(tc_dct_small, test_dct_EMx_new_small);
  suite_add_tcase(s, tc_dct_small);


  tc_dct_large = tcase_create("dct_large");
  tcase_add_checked_fixture(tc_dct_large, setup_large, teardown_large);
  tcase_add_test(tc_dct_large,test_dct_MtEt_EMx_large);
  tcase_add_test(tc_dct_large, test_dct_MtEty_large);
  tcase_add_test(tc_dct_large, test_dct_EMx_large);

  suite_add_tcase(s, tc_dct_large);

  tc_dct_pure = tcase_create("dct_pure");
  tcase_add_checked_fixture(tc_dct_pure, setup_pure_dct, teardown_pure);
  tcase_add_test(tc_dct_pure, test_dct_MtEt_EMx_pure);
  tcase_add_test(tc_dct_pure, test_dct_EMx_new_pure);
  tcase_add_test(tc_dct_pure, test_dct_MtEty_pure);

  suite_add_tcase(s, tc_dct_pure);




  return s;

}
