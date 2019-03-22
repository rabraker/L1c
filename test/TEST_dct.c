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

#include "dct.h"

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
DctData *dctd;
// _small_data;
// DctData *dct_large_data;
// DctData *dctd;



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

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  dct_dat->MtEty_EMx_act = malloc_double(dct_dat->Nx);
  dct_dat->MtEty_act = malloc_double(dct_dat->Nx);
  dct_dat->EMx_act = malloc_double(dct_dat->Nx);

  dct_setup(dct_dat->Nx, dct_dat->Npix, dct_dat->pix_idx);

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
  dct_destroy();


#ifdef _USEMKL_
  mkl_free_buffers();
#endif

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
  dct_MtEt_EMx(dctd->x_in, dctd->MtEty_EMx_act);

  ck_assert_double_array_eq_tol(dctd->Nx, dctd->MtEty_EMx_exp,
                                dctd->MtEty_EMx_act, TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST(test_dct_MtEty)
{
  dct_MtEty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->Nx, dctd->MtEty_exp,
                                dctd->MtEty_act, TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST(test_dct_EMx)
{
  dct_EMx(dctd->x_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->Ny, dctd->EMx_exp,
                                dctd->EMx_act, TOL_DOUBLE_SUPER);

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
  suite_add_tcase(s, tc_dct_small);

  tc_dct_large = tcase_create("dct_large");
  tcase_add_checked_fixture(tc_dct_large, setup_large, teardown_large);
  tcase_add_test(tc_dct_large,test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_large, test_dct_MtEty);
  tcase_add_test(tc_dct_large, test_dct_EMx);
  suite_add_tcase(s, tc_dct_large);

  tc_dct_pure = tcase_create("dct_pure");
  tcase_add_checked_fixture(tc_dct_pure, setup_pure, teardown_pure);
  tcase_add_test(tc_dct_pure, test_dct_MtEt_EMx);
  tcase_add_test(tc_dct_pure, test_dct_EMx);
  tcase_add_test(tc_dct_pure, test_dct_MtEty);
  suite_add_tcase(s, tc_dct_pure);

  return s;

}
