#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#include "l1c_common.h"

#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"
#include "l1c_math.h"


START_TEST(test_l1c_daxpy_z)
{
  l1c_int N = 6;
  double a = 3;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double y[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  double z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  double z_exp[] = {5.0, 8.0, 11.0, 14.0, 17.0, 20.0};


  l1c_daxpy_z(N, a, x, y, z);

  ck_assert_double_array_eq_tol(N, z_exp, z, TOL_DOUBLE);
}
END_TEST

START_TEST(test_l1c_dxmuly_z)
{
  l1c_int N = 6;

  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double y[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  double z[] = {0.0, 0.0, 0.0, 0.0, 0.0, .0};
  double z_exp[] = {2.0, 4.0, 6.0, 8.0, 10.0, 12.0};
  l1c_dxmuly_z(N, x, y, z);

  ck_assert_double_array_eq_tol(N, z_exp, z, TOL_DOUBLE);
}
END_TEST




START_TEST(test_l1c_dsum)
{

  l1c_int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double sum_exp = 21.0;
  double sum_x = l1c_dsum(N, x);

  ck_assert_double_eq_tol(sum_exp, sum_x, TOL_DOUBLE);
}
END_TEST

START_TEST(test_l1c_dnorm1)
{

  l1c_int N = 6;
  double x[] = {1.0, 2.0, -3.0, -4.0, 5.0, 6.0};
  double sum_exp = 21.0;
  double sum_x = l1c_dnorm1(N, x);

  ck_assert_double_eq_tol(sum_exp, sum_x, TOL_DOUBLE);
}
END_TEST

START_TEST(test_l1c_dlogsum)
{
  l1c_int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double alpha = 3.0;
  double logsum_x_exp = 13.170924944018758;
  double logsum_x = l1c_dlogsum(N, alpha, x);

  ck_assert_double_eq_tol(logsum_x_exp, logsum_x, TOL_DOUBLE);
  // setup_vectors_SC();
}
END_TEST




Suite *l1c_math_suite(void)
{
  Suite *s;

  TCase *tc_l1c_math;

  s = suite_create("l1c_math");


  tc_l1c_math = tcase_create("l1qc_math_funs");

  tcase_add_test(tc_l1c_math, test_l1c_dxmuly_z);
  tcase_add_test(tc_l1c_math, test_l1c_daxpy_z);
  tcase_add_test(tc_l1c_math, test_l1c_dsum);
  tcase_add_test(tc_l1c_math, test_l1c_dnorm1);
  tcase_add_test(tc_l1c_math, test_l1c_dlogsum);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_l1c_math);

  return s;

}
