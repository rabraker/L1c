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

// #include "test_data.h"
#include "cgsolve.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "cJSON.h"
#include "json_utils.h"

#include "l1qc_common.h"
#include "bregman.h"
#include "check_utils.h"

cJSON *test_data_json;




START_TEST(test_breg_diffX){
  BregFuncs bfuncs = breg_get_functions();

  char fpath[] = "test_data/bregman_diff.json";

  double *A=NULL, *Dx=NULL, *Dx_exp=NULL;

  l1c_int N=0, m=0, n=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "A", &A, &N);
  status +=extract_json_double_array(test_data_json, "dx", &Dx_exp, &N);
  status +=extract_json_int(test_data_json, "n", &n);
  status +=extract_json_int(test_data_json, "m", &m);

  Dx = malloc_double(N);
  if ( (!Dx) | status){
    perror("error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    Dx[i] =0;
  }

  bfuncs.breg_Dx(n, m, A, Dx);

  ck_assert_double_array_eq_tol(N, Dx_exp, Dx, TOL_DOUBLE_SUPER*10);


  free_double(A);
  free_double(Dx);
  free_double(Dx_exp);


  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_breg_diffY){

  BregFuncs bfuncs = breg_get_functions();
  char fpath[] = "test_data/bregman_diff.json";

  double *A=NULL, *Dy=NULL, *Dy_exp=NULL;

  l1c_int N=0, m=0, n=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }


  status +=extract_json_double_array(test_data_json, "A", &A, &N);
  status +=extract_json_double_array(test_data_json, "dy", &Dy_exp, &N);
  status +=extract_json_int(test_data_json, "n", &n);
  status +=extract_json_int(test_data_json, "m", &m);

  Dy = malloc_double(N);
  if ( (!Dy) | status){
    perror("error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    Dy[i] =0;
  }
  bfuncs.breg_Dy(n, m, A, Dy);

  ck_assert_double_array_eq_tol(N, Dy_exp, Dy, TOL_DOUBLE_SUPER*10);


  free_double(A);
  free_double(Dy);
  free_double(Dy_exp);


  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_breg_shrink1){

  BregFuncs bfuncs = breg_get_functions();
  char fpath[] = "test_data/bregman.json";

  double *x=NULL, *x_shrunk=NULL, *x_shrunk_exp=NULL;
  double gamma;
  l1c_int N=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "x_shrunk", &x_shrunk_exp, &N);
  status +=extract_json_double(test_data_json, "gamma", &gamma);

  x_shrunk = malloc_double(N);
  if ( (!x_shrunk) | status){
    perror("error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    x_shrunk[i] =0;
  }
  bfuncs.breg_shrink1(N, x, x_shrunk, gamma);

  ck_assert_double_array_eq_tol(N, x_shrunk_exp, x_shrunk, TOL_DOUBLE_SUPER*10);


  free_double(x);
  free_double(x_shrunk);
  free_double(x_shrunk_exp);


  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_breg_mxpy_z){

  BregFuncs bfuncs = breg_get_functions();
  char fpath[] = "test_data/bregman.json";

  double *x=NULL, *y=NULL, *z_exp=NULL, *z=NULL;
  double gamma;
  l1c_int N=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_get_gradient\n");
    ck_abort();
  }

  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "y", &y, &N);
  status +=extract_json_double_array(test_data_json, "z", &z_exp, &N);
  status +=extract_json_double(test_data_json, "gamma", &gamma);

  z = malloc_double(N);
  if ( (!z) | status){
    perror("error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    z[i] =1;
  }
  bfuncs.breg_mxpy_z(N, x, y, z);

  ck_assert_double_array_eq_tol(N, z_exp, z, TOL_DOUBLE_SUPER*10);


  free_double(x);
  free_double(y);
  free_double(z);
  free_double(z_exp);


  free_json_text_data();
  cJSON_Delete(test_data_json);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST

/* Add all the test cases to our suite
 */
Suite *bregman_suite(void)
{
  Suite *s;

  TCase *tc_diff;
  s = suite_create("bregman");
  tc_diff = tcase_create("diffX");

  tcase_add_test(tc_diff, test_breg_diffX);
  tcase_add_test(tc_diff, test_breg_diffY);
  tcase_add_test(tc_diff, test_breg_shrink1);
  tcase_add_test(tc_diff, test_breg_mxpy_z);

  suite_add_tcase(s, tc_diff);

  return s;

}
