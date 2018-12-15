/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */

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

/* Global variables which hold data contained in
   test_data_ss_ff.h
*/

// double *y_vec_tmp, *DC_y_tmp;

static int load_EMx_data(int *Nx0, double **x0, int *Nx1, double **x1, int *Nidx,
                         int **pix_idx, char *fpath){

  cJSON *test_data_json;

  load_file_to_json(fpath, &test_data_json);


  if (extract_json_double_array(test_data_json, "x0", x0, Nx0) ){
    perror("Error Loading x\n");
    return 1;
  }
  if (extract_json_double_array(test_data_json, "x1", x1, Nx1) ){
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

 // end:
 //  return 1;
}



START_TEST(test_dct_MtEt_EMtx_small)
{
  /* Test the multiplication (EM)^T * (E*M) * x */
  char fpath[] = "test_data/dct_small_MtEt_EMx.json";
  int *pix_idx;
  double *x_exp, *x0, *x_act;
  int Nx0, Nx1, Npix = 0;

  if (load_EMx_data(&Nx0, &x0, &Nx1, &x_exp, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data\n");
  }

  printf("Nx0 = %d, Nx1=%d\n", Nx0, Nx1);
  dct_setup(Nx0, Npix, pix_idx);

  x_act = dct_MtEt_EMx_new(x0);

  ck_assert_double_array_eq_tol(Nx0, x_exp, x_act, TOL_DOUBLE);

  free(x_exp);
  dct_destroy();
}
END_TEST


START_TEST(test_dct_MtEty_small)
{
  char fpath[] = "test_data/dct_small_MtEty.json";
  int *pix_idx;
  double *x_exp, *y, *x_act;
  int Nx, Ny, Npix = 0;

  if (load_EMx_data(&Nx, &x_exp, &Ny, &y, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data\n");
  }

  printf("Nx = %d, Ny=%d\n", Nx, Ny);
  dct_setup(Nx, Ny, pix_idx);

  x_act = dct_MtEty(y);

  ck_assert_double_array_eq_tol(Nx, x_exp, x_act, TOL_DOUBLE);

  free(y);
  dct_destroy();
}
END_TEST


START_TEST(test_dct_EMx_small)
{
  char fpath[] = "test_data/dct_small_EMx.json";

  int *pix_idx;

  double *x, *y_exp, *y_act;
  int Nx, Ny, Npix = 0;
  if (load_EMx_data(&Nx, &x, &Ny, &y_exp, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data\n");
  }

  printf("Nx = %d, Ny=%d\n", Nx, Ny);
  dct_setup(Nx, Ny, pix_idx);
  dct_load_x(x);

  y_act = dct_EMx();

  ck_assert_double_array_eq_tol(Ny, y_exp, y_act, TOL_DOUBLE);

  free(y_exp);
  dct_destroy();
}
END_TEST

START_TEST(test_dct_EMx_new_small)
{
  char fpath[] = "test_data/dct_small_EMx.json";

  int *pix_idx;

  double *x, *x_new, *y_exp, *y_act;
  int Nx, Ny, Npix, i = 0;
  if (load_EMx_data(&Nx, &x, &Ny, &y_exp, &Npix, &pix_idx, fpath)){
    ck_abort_msg("Errory Loading test data\n");
  }
  // x from json loader is not properly aligned.
  // load the data into x_new, which has been properly allocated.
  x_new = fftw_alloc_real(Nx);
  for (i=0; i<Nx; i++){
    x_new[i] = x[i];
  }

  printf("Nx = %d, Ny=%d\n", Nx, Ny);
  dct_setup(Nx, Ny, pix_idx);
  //dct_load_x(x);

  y_act = dct_EMx_new(x_new);

  ck_assert_double_array_eq_tol(Ny, y_exp, y_act, TOL_DOUBLE);

  free(y_exp);
  dct_destroy();
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

  tcase_add_test(tc_core, test_dct_EMx_small);
  tcase_add_test(tc_core, test_dct_EMx_new_small);

  tcase_add_test(tc_core, test_dct_MtEt_EMtx_small);
  tcase_add_test(tc_core, test_dct_MtEty_small);



  suite_add_tcase(s, tc_core);

  return s;

}
