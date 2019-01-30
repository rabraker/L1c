/*
This is a test suite for the fourier integration functions contained in
../src/ss_fourier_functions.c. The data used to test these functions is
contained in the header file test_data_ss_ff.h, which defines several global
variables. The header file is generated from the matlab script called
generate_test_data.m

This test suite uses the libcheck framework. On my computer, this got installed
into /usr/local/lib, which was not by default found by the system. Thus, I have
to do
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib


  https://libcheck.github.io/check/doc/check_html/check_4.html#No-Fork-Mode

 */


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#include "l1c_common.h"

#include <cjson/cJSON.h>
#include "json_utils.h"
#include "vcl_math.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"



START_TEST(test_vcl_sum)
{
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double *x = malloc_double(12);
  for (int i =0; i<12; i++){
    x[i] = x_[i];
  }

  double sum_exp0 = 6.0;
  double sum_exp1 = 36.0;
  double sum_exp2 = 45.0;
  double sum_exp3 = 55.0;
  double sum_exp4 = 66.0;
  double sum_exp5 = 78.0;
  double sum_x = 0.0;

  sum_x = vcl_sum(3, x);
  ck_assert_double_eq_tol(sum_exp0, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(8, x);
  ck_assert_double_eq_tol(sum_exp1, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(9, x);
  ck_assert_double_eq_tol(sum_exp2, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(10, x);
  ck_assert_double_eq_tol(sum_exp3, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(11, x);
  ck_assert_double_eq_tol(sum_exp4, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(12, x);
  ck_assert_double_eq_tol(sum_exp5, sum_x, TOL_DOUBLE);

  free_double(x);
}
END_TEST


START_TEST(test_vcl_logsum)
{
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double *x = malloc_double(12);
  for (int i =0; i<12; i++){
    x[i] = x_[i];
  }
  double logsum_exp0 = 5.087596335232384;
  double logsum_exp1 = 19.393501212090126;
  double logsum_exp2 = 22.689338078094455;
  double logsum_exp3 = 26.090535459756609;
  double logsum_exp4 = 29.587043021223089;
  double logsum_exp5 = 33.170561959679198;

  double logsum_x = 0.0;

  double alpha = 3.0;

  logsum_x = vcl_logsum(3, alpha, x);
  ck_assert_double_eq_tol(logsum_exp0, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(8, alpha, x);
  ck_assert_double_eq_tol(logsum_exp1, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(9, alpha, x);
  ck_assert_double_eq_tol(logsum_exp2, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(10, alpha, x);
  ck_assert_double_eq_tol(logsum_exp3, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(11, alpha, x);
  ck_assert_double_eq_tol(logsum_exp4, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(12, alpha, x);
  ck_assert_double_eq_tol(logsum_exp5, logsum_x, TOL_DOUBLE);

  free_double(x);
}
END_TEST




Suite *vcl_math_suite(void)
{
  Suite *s;

  TCase  *tc_mathfuns;
  s = suite_create("vcl_suite");


  tc_mathfuns = tcase_create("vcl_math");
  tcase_add_test(tc_mathfuns, test_vcl_sum);
  tcase_add_test(tc_mathfuns, test_vcl_logsum);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_mathfuns);

  return s;

}
