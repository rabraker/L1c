/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */

#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#ifdef _USEMKL_  //Disable the entire thing if no mkl. See end of file.


// #include "test_data.h"
#include "dct_mkl.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "cJSON.h"
#include "json_utils.h"
#include "check_utils.h"
#include "l1c_common.h"

#define TOL_LARGE_DCT 1e-10


/* Global variables which hold data contained in
   test_data_ss_ff.h
*/

cJSON *test_data_json;


START_TEST(test_dctmkl_MtEt_EMx_small_rand)
{
  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *MtEt_EMx_exp, *x_in, *MtEt_EMx_act;
  l1c_int Nx, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "MtEt_EMx", &MtEt_EMx_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  dctmkl_setup(Nx, Npix, (MKL_INT*)pix_idx);

  MtEt_EMx_act = malloc_double(Nx);

  dctmkl_MtEt_EMx_new(x_in, MtEt_EMx_act);

  ck_assert_double_array_eq_tol(Nx, MtEt_EMx_exp, MtEt_EMx_act, TOL_DOUBLE_SUPER);

  free_double(x_in);
  free_double(MtEt_EMx_act);
  free_double(MtEt_EMx_exp);
  free(pix_idx);
  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();
}
END_TEST


START_TEST(test_dctmkl_MtEty_small_rand)
{

  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *MtEty_exp, *y_in, *MtEty_act;
  l1c_int Nx, Ny, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "y_in", &y_in, &Ny);
  status +=extract_json_double_array(test_data_json, "MtEty", &MtEty_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  MtEty_act = malloc_double(Nx);
  ck_assert_int_eq(Ny, Npix);

  dctmkl_setup(Nx, Ny, pix_idx);

  dctmkl_MtEty(y_in, MtEty_act);

  ck_assert_double_array_eq_tol(Nx, MtEty_exp, MtEty_act, TOL_DOUBLE_SUPER);

  free_double(y_in);
  free_double(MtEty_act);
  free_double(MtEty_exp);
  free(pix_idx);
  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();
}
END_TEST

START_TEST(test_dctmkl_EMx_new_small_rand)
{

  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *EMx_exp, *x_in, *x_in_aligned, *EMx_act;
  l1c_int Nx, Ny, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "EMx", &EMx_exp, &Ny);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(Ny, Npix);

  dctmkl_setup(Nx, Ny, pix_idx);

  x_in_aligned = malloc_double(Nx);
  EMx_act = malloc_double(Nx);
  for (l1c_int i=0; i<Nx; i++){
    x_in_aligned[i] = x_in[i];
  }

  dctmkl_EMx_new(x_in_aligned, EMx_act);

  ck_assert_double_array_eq_tol(Ny, EMx_exp, EMx_act, TOL_DOUBLE_SUPER);

  free_double(x_in);
  free_double(EMx_act);
  free_double(EMx_exp);
  free(pix_idx);
  free_double(x_in_aligned);
  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();
}
END_TEST



START_TEST(test_dctmkl_MtEt_EMx_large)
{
  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *MtEt_EMx_exp, *x_in, *MtEt_EMx_act;
  l1c_int Nx, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "MtEt_EMx", &MtEt_EMx_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  dctmkl_setup(Nx, Npix, pix_idx);
  MtEt_EMx_act = malloc_double(Nx);

  dctmkl_MtEt_EMx_new(x_in, MtEt_EMx_act);

  ck_assert_double_array_eq_tol(Nx, MtEt_EMx_exp, MtEt_EMx_act,  TOL_LARGE_DCT);

  free_double(x_in);
  free_double(MtEt_EMx_exp);
  free(pix_idx);
  free_double(MtEt_EMx_act);

  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();

}
END_TEST


START_TEST(test_dctmkl_MtEty_large)
{

  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *MtEty_exp, *y_in, *MtEty_act;
  l1c_int Nx, Ny, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "y_in", &y_in, &Ny);
  status +=extract_json_double_array(test_data_json, "MtEty", &MtEty_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  ck_assert_int_eq(Ny, Npix);

  dctmkl_setup(Nx, Ny, pix_idx);
  MtEty_act = malloc_double(Nx);
  dctmkl_MtEty(y_in, MtEty_act);

  ck_assert_double_array_eq_tol(Nx, MtEty_exp, MtEty_act, TOL_LARGE_DCT);

  free_double(y_in);
  free_double(MtEty_exp);
  free(pix_idx);
  free_double(MtEty_act);
  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();

}
END_TEST


START_TEST(test_dctmkl_EMx_large)
{

  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *EMx_exp, *x_in, *EMx_act;
  l1c_int Nx, Ny, Npix, status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "EMx", &EMx_exp, &Ny);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(Ny, Npix);

  dctmkl_setup(Nx, Ny, pix_idx);
  EMx_act = malloc_double(Ny);

  dctmkl_EMx_new(x_in, EMx_act);

  ck_assert_double_array_eq_tol(Ny, EMx_exp, EMx_act, TOL_LARGE_DCT);

  free_double(x_in);
  free_double(EMx_exp);
  free_double(EMx_act);
  free(pix_idx);

  dctmkl_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
  mkl_free_buffers();
}
END_TEST




/* Add all the test cases to our suite
 */
Suite *dct_mkl_suite(void)
{
  Suite *s;
  TCase *tc_small, *tc_large;
  s = suite_create("dct_mkl");
  tc_large = tcase_create("dct_mkl_large");
  tc_small = tcase_create("dct_mkl_small");

  tcase_add_test(tc_small, test_dctmkl_MtEt_EMx_small_rand);
  tcase_add_test(tc_small, test_dctmkl_EMx_new_small_rand);
  tcase_add_test(tc_small, test_dctmkl_MtEty_small_rand);

  tcase_add_test(tc_large, test_dctmkl_MtEt_EMx_large);
  tcase_add_test(tc_large, test_dctmkl_EMx_large);
  tcase_add_test(tc_large, test_dctmkl_MtEty_large);

  suite_add_tcase(s, tc_small);
  suite_add_tcase(s, tc_large);

  return s;

}

#else //Just a dummy for when mkl is disabled.
Suite *dct_mkl_suite(void)
{
  Suite *s;
  // TCase *tc_core;
  s = suite_create("dct_mkl");

  return s;
}

#endif
