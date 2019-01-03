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
#include <fftw3.h>

// #include "test_data.h"
#include "dct.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "cJSON.h"
#include "json_utils.h"
#include "check_utils.h"
#include "l1qc_common.h"

#ifdef _USEMKL_
#define TOL_LARGE_DCT 1e-10
#else
#define TOL_LARGE_DCT TOL_DOUBLE_SUPER
#endif

/* Global variables which hold data contained in
   test_data_ss_ff.h
*/

cJSON *test_data_json;


START_TEST(test_dct_MtEt_EMx_small_rand)
{
  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *MtEt_EMx_exp, *x_in, *MtEt_EMx_act;
  l1c_int Nx, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "MtEt_EMx", &MtEt_EMx_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  dct_setup(Nx, Npix, pix_idx);

  MtEt_EMx_act = fftw_alloc_real(Nx);
  dct_MtEt_EMx_new(x_in, MtEt_EMx_act);

  ck_assert_double_array_eq_tol(Nx, MtEt_EMx_exp, MtEt_EMx_act, TOL_DOUBLE_SUPER);

  fftw_free(x_in);
  fftw_free(MtEt_EMx_act);
  free_double(MtEt_EMx_exp);
  free(pix_idx);
  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST


START_TEST(test_dct_MtEty_small_rand)
{

  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *MtEty_exp, *y_in, *MtEty_act;
  l1c_int Nx, Ny, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "y_in", &y_in, &Ny);
  status +=extract_json_double_array(test_data_json, "MtEty", &MtEty_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  MtEty_act = fftw_alloc_real(Nx);
  ck_assert_int_eq(Ny, Npix);

  dct_setup(Nx, Ny, pix_idx);

  dct_MtEty(y_in, MtEty_act);

  ck_assert_double_array_eq_tol(Nx, MtEty_exp, MtEty_act, TOL_DOUBLE_SUPER);

  fftw_free(y_in);
  fftw_free(MtEty_act);
  free_double(MtEty_exp);
  free(pix_idx);
  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST

START_TEST(test_dct_EMx_new_small_rand)
{

  char fpath[] = "test_data/dct_small_rand.json";
  l1c_int *pix_idx;
  double *EMx_exp, *x_in, *x_in_aligned, *EMx_act;
  l1c_int Nx, Ny, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "EMx", &EMx_exp, &Ny);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(Ny, Npix);

  dct_setup(Nx, Ny, pix_idx);

  x_in_aligned = fftw_alloc_real(Nx);
  EMx_act = fftw_alloc_real(Nx);
  for (l1c_int i=0; i<Nx; i++){
    x_in_aligned[i] = x_in[i];
  }

  dct_EMx_new(x_in_aligned, EMx_act);

  ck_assert_double_array_eq_tol(Ny, EMx_exp, EMx_act, TOL_DOUBLE_SUPER);

  fftw_free(x_in);
  fftw_free(EMx_act);
  free_double(EMx_exp);
  free(pix_idx);
  fftw_free(x_in_aligned);
  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST



/* ---------------------------------------------------- */
START_TEST(test_dct_MtEt_EMx_large)
{
  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *MtEt_EMx_exp, *x_in, *MtEt_EMx_act;
  l1c_int Nx, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "MtEt_EMx", &MtEt_EMx_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  dct_setup(Nx, Npix, pix_idx);
  MtEt_EMx_act = malloc_double(Nx);

  dct_MtEt_EMx_new(x_in, MtEt_EMx_act);

  ck_assert_double_array_eq_tol(Nx, MtEt_EMx_exp, MtEt_EMx_act,  TOL_LARGE_DCT);

  fftw_free(x_in);
  free_double(MtEt_EMx_exp);
  free(pix_idx);
  free_double(MtEt_EMx_act);

  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST


START_TEST(test_dct_MtEty_large)
{

  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *MtEty_exp, *y_in, *MtEty_act;
  l1c_int Nx, Ny, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "y_in", &y_in, &Ny);
  status +=extract_json_double_array(test_data_json, "MtEty", &MtEty_exp, &Nx);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }

  ck_assert_int_eq(Ny, Npix);

  dct_setup(Nx, Ny, pix_idx);
  MtEty_act = malloc_double(Nx);
  dct_MtEty(y_in, MtEty_act);

  ck_assert_double_array_eq_tol(Nx, MtEty_exp, MtEty_act, TOL_LARGE_DCT);

  fftw_free(y_in);
  free_double(MtEty_exp);
  free(pix_idx);
  free_double(MtEty_act);
  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST

START_TEST(test_dct_EMx_large)
{

  char fpath[] = "test_data/dct_large.json";
  l1c_int *pix_idx;
  double *EMx_exp, *x_in, *EMx_act;
  l1c_int Nx, Ny, Npix;
  int status = 0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_dct_MtEt_large\n");
    ck_abort();
  }

  status +=extract_json_double_array_fftw(test_data_json, "x_in", &x_in, &Nx);
  status +=extract_json_double_array(test_data_json, "EMx", &EMx_exp, &Ny);
  status +=extract_json_int_array(test_data_json, "pix_idx", &pix_idx, &Npix);

  if (status){
    perror("Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(Ny, Npix);

  dct_setup(Nx, Ny, pix_idx);
  EMx_act = malloc_double(Ny);

  dct_EMx_new(x_in, EMx_act);

  ck_assert_double_array_eq_tol(Ny, EMx_exp, EMx_act, TOL_LARGE_DCT);

  fftw_free(x_in);
  free_double(EMx_exp);
  free_double(EMx_act);
  free(pix_idx);

  dct_destroy(); //will free MtEty_act.

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST



static int load_EMx_data(l1c_int *Nx0, double **x0, l1c_int *Nx1, double **x1, l1c_int *Nidx,
                         l1c_int **pix_idx, char *fpath){

  if (load_file_to_json(fpath, &test_data_json) ){
    return 1;
  }
  if (extract_json_double_array_fftw(test_data_json, "x0", x0, Nx0) ){
    perror("Error Loading x\n");
    return 1;
  }
  if (extract_json_double_array_fftw(test_data_json, "x1", x1, Nx1) ){
    perror("Error Loading y_exp\n");
    return 1;
  }
  if (extract_json_int_array(test_data_json, "pix_idx", pix_idx, Nidx) ){
    perror("Error Loading pix_idx \n");
    goto end1;
  }

  /* Sanity check */
  if ( (*Nx1 != *Nidx) && (*Nx0 != *Nx1)){
    perror("Error: Array size mismatch. Aborting\n");
    goto end2; // We allocated all, but their sizes don't match.
  }

  return 0;

 end2:
  free(*pix_idx);
  goto end1;
 end1:
  free(*x1);
  goto end0;
 end0:
  free(*x0);
  return 1;

}


START_TEST(test_dct_MtEt_EMx_small)
{
  /* Test the multiplication (EM)^T * (E*M) * x */
  char fpath[] = "test_data/dct_small_MtEt_EMx.json";
  l1c_int *pix_idx;
  double *x_exp, *x0, *x_act;
  l1c_int Nx0, Nx1, Npix = 0;

  if (load_EMx_data(&Nx0, &x0, &Nx1, &x_exp, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data (in test_dct_MtEt_EMx\n");
  }

  dct_setup(Nx0, Npix, pix_idx);
  x_act = malloc_double(Nx0);

  dct_MtEt_EMx_new(x0, x_act);

  ck_assert_double_array_eq_tol(Nx0, x_exp, x_act, TOL_DOUBLE_SUPER);

  fftw_free(x_exp);
  fftw_free(x0);
  free(pix_idx);
  dct_destroy();
  free_double(x_act);
  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST


START_TEST(test_dct_MtEty_small)
{
  char fpath[] = "test_data/dct_small_MtEty.json";
  l1c_int *pix_idx;
  double *x_exp, *y, *x_act;
  l1c_int Nx, Ny, Npix = 0;

  if (load_EMx_data(&Nx, &x_exp, &Ny, &y, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data (in test_dct_MtEty) \n");
  }

  dct_setup(Nx, Ny, pix_idx);
  x_act = malloc_double(Nx);
  dct_MtEty(y, x_act);

  ck_assert_double_array_eq_tol(Nx, x_exp, x_act, TOL_DOUBLE_SUPER);

  fftw_free(y);
  fftw_free(x_exp);
  dct_destroy();
  free_double(x_act);

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST




START_TEST(test_dct_EMx_new_small)
{
  char fpath[] = "test_data/dct_small_EMx.json";

  l1c_int *pix_idx;

  double *x, *x_new, *y_exp, *y_act;
  l1c_int Nx, Ny, Npix, i = 0;
  if (load_EMx_data(&Nx, &x, &Ny, &y_exp, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data\n");
  }
  // x from json loader is not properly aligned.
  // load the data into x_new, which has been properly allocated.
  x_new = fftw_alloc_real(Nx);
  for (i=0; i<Nx; i++){
    x_new[i] = x[i];
  }

  dct_setup(Nx, Ny, pix_idx);
  y_act = malloc_double(Ny);

  dct_EMx_new(x_new, y_act);

  ck_assert_double_array_eq_tol(Ny, y_exp, y_act, TOL_DOUBLE_SUPER);

  fftw_free(y_exp);
  fftw_free(x_new);
  fftw_free(x);
  free(pix_idx);
  dct_destroy();
  free_double(y_act);

  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif
}
END_TEST




/* Add all the test cases to our suite
 */
Suite *dct_suite(void)
{
  Suite *s;

  TCase *tc_core;
  s = suite_create("dct");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core,test_dct_MtEt_EMx_small_rand);
  tcase_add_test(tc_core, test_dct_EMx_new_small_rand);
  tcase_add_test(tc_core, test_dct_MtEty_small_rand);

  tcase_add_test(tc_core,test_dct_MtEt_EMx_large);
  tcase_add_test(tc_core, test_dct_EMx_large);
  tcase_add_test(tc_core, test_dct_MtEty_large);

  tcase_add_test(tc_core, test_dct_EMx_new_small);
  tcase_add_test(tc_core, test_dct_MtEt_EMx_small);
  tcase_add_test(tc_core, test_dct_MtEty_small);



  suite_add_tcase(s, tc_core);

  return s;

}
