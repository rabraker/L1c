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

// #include "test_data.h"
#include "cgsolve.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include <cjson/cJSON.h>
#include "json_utils.h"

#include "l1c_common.h"
#include "bregman.h"
#include "check_utils.h"

cJSON *test_data_json;


START_TEST(test_breg_shrink1){

  BregFuncs bfuncs = breg_get_functions();
  char fpath[] = "test_data/bregman.json";

  double *x=NULL, *x_shrunk=NULL, *x_shrunk_exp=NULL;
  double gamma;
  l1c_int N=0, status=0;

  if (load_file_to_json(fpath, &test_data_json)){
    perror("Error loading data in test_breg_shrink1\n");
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
    perror("Error loading data in test_breg_mxpy_z\n");
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

  tcase_add_test(tc_diff, test_breg_shrink1);
  tcase_add_test(tc_diff, test_breg_mxpy_z);

  suite_add_tcase(s, tc_diff);

  return s;

}
