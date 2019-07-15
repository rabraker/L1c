#include "config.h"

#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#include "l1c.h"

#include <cjson/cJSON.h>
#include "json_utils.h"
#include "l1qc_newton.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"
#include "l1c_math.h"



START_TEST(test_l1c_init_vec)
{
  l1c_int N = 6;
  double a = 3.5;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  l1c_init_vec(N, x, a);

  for(int i=0; i<N; i++){
    ck_assert_double_eq_tol(a, x[i], TOL_DOUBLE_SUPER);
  }

}
END_TEST

START_TEST(test_l1c_daxpy_z2)
{
  l1c_int N = 256*256;

  double *xx = l1c_malloc_double(2*N);
  double *yy = l1c_malloc_double(2*N);
  double *zz_exp = l1c_malloc_double(2*N);
  double *zz = l1c_malloc_double(2*N);

  double alp = 4.5918;
  for (int i=0; i<2*N; i++){
    xx[i] = ((double) rand()) / 256 / 256;
    yy[i] = ((double) rand()) / 256 / 256;
  }

  for (int i=0; i<2*N; i++){
    zz_exp[i] = alp * xx[i] + yy[i];
  }

  l1c_daxpy_z(N, alp, xx, yy, zz);
  l1c_daxpy_z(N, alp, xx+N, yy+N, zz+N);

  ck_assert_double_array_eq_tol(N, zz_exp, zz, TOL_DOUBLE_SUPER);

  l1c_free_double(xx);
  l1c_free_double(yy);
  l1c_free_double(zz_exp);
  l1c_free_double(zz);

}
END_TEST


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

START_TEST(test_l1c_daxpby_z)
{
  l1c_int N = 6;
  double a = 3;
  double b = 2;

  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double y[] = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
  double z[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  double z_exp[] = {7.0, 10.0, 13.0, 16.0, 19.0, 22.0};


  l1c_daxpby_z(N, a, x, b, y, z);

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

START_TEST(test_l1c_dnrm2_err)
{
  l1c_int N = 6;
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double y_[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  double *x = l1c_malloc_double(N);
  double *y = l1c_malloc_double(N);
  if (!x || !y){
    ck_abort_msg("Out of memory\n");
  }
  for (int i=0; i<N; i++){
    x[i] = x_[i];
    y[i] = y_[i];
  }

  double nrm_x_y = l1c_dnrm2_err(N, x, y);

  ck_assert_double_eq_tol(nrm_x_y, sqrt(N), TOL_DOUBLE_SUPER);

  l1c_free_double(x);
  l1c_free_double(y);
}
END_TEST

START_TEST(test_l1c_dnrm2_rel_err)
{
  l1c_int N = 6;
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double y_[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  double *x = l1c_malloc_double(N);
  double *y = l1c_malloc_double(N);
  if (!x || !y){
    ck_abort_msg("Out of memory\n");
  }
  for (int i=0; i<N; i++){
    x[i] = x_[i];
    y[i] = y_[i];
  }

  double nrm_x_y = l1c_dnrm2_rel_err(N, x, y);

  ck_assert_double_eq_tol(nrm_x_y, sqrt(N)/sqrt(55.0), TOL_DOUBLE_SUPER);

  l1c_free_double(x);
  l1c_free_double(y);
}
END_TEST

START_TEST(test_l1c_max_vec)
{
  l1c_int N = 6;
  double xmax;
  double x_[] = {1.0, -2.0, -3.0, 4.0, 5.0, 6.0};

  double *x=l1c_malloc_double(N);
  if (!x){
    ck_abort_msg("Out of memory\n");
  }
  for (int i=0; i<N; i++){
    x[i] = x_[i];
  }

  xmax = l1c_max_vec(0, x);
  //ck_assert_double_nan(xmax);
  ck_assert_msg(isnan(xmax),  "Assertion is NaN failed: xmax == %f\n", xmax);


  xmax = l1c_max_vec(1, x);
  ck_assert_double_eq(1.0, xmax);

  xmax = l1c_max_vec(N, x);
  ck_assert_double_eq(6.0, xmax);

  l1c_free_double(x);
}
END_TEST


START_TEST(test_l1c_abs_vec)
{
  l1c_int N = 6;
  double x_[]    = {1.0, -2.0, -3.0, 4.0, 5.0, 6.0};
  double x_exp[] = {1.0,  2.0,  3.0, 4.0, 5.0, 6.0};

  double *x = l1c_malloc_double(N);
  double *xabs = l1c_malloc_double(N);
  if (!x || !xabs){
    ck_abort_msg("Out of memory\n");
  }
  for (int i=0; i<N; i++){
    x[i] = x_[i];
    xabs[i] = 10;
  }

  l1c_abs_vec(N, x, xabs);

  ck_assert_double_array_eq_tol(N, xabs, x_exp, TOL_DOUBLE_SUPER);

  l1c_free_double(x);
  l1c_free_double(xabs);
}
END_TEST

START_TEST(l1c_math_max)
{
  double a=5, b=6, mx = 0;

  mx = max(a, b);
  ck_assert_double_eq(mx, b);

  mx = max(a, -b);
  ck_assert_double_eq(mx, a);
}
END_TEST

START_TEST(l1c_math_min)
{
  double a=5, b=6, mx = 0;

  mx = min(a, b);
  ck_assert_double_eq(mx, a);

  mx = min(a, -b);
  ck_assert_double_eq(mx, -b);

}
END_TEST

START_TEST(l1c_math_imax)
{
  l1c_int a=5, b=6, mx = 0;

  mx = imax(a, b);
  ck_assert_int_eq(mx, b);

  mx = imax(a, -b);
  ck_assert_int_eq(mx, a);
}
END_TEST

START_TEST (l1c_math_imin)
{
  l1c_int a=5, b=6, mx = 0;

  mx = imin(a, b);
  ck_assert_double_eq(mx, a);

  mx = imin(a, -b);
  ck_assert_double_eq(mx, -b);

}
END_TEST


Suite *l1c_math_suite(void)
{
  Suite *s;

  TCase *tc_l1c_math;

  s = suite_create("l1c_math");


  tc_l1c_math = tcase_create("l1qc_math_funs");

  tcase_add_test(tc_l1c_math, test_l1c_init_vec);
  tcase_add_test(tc_l1c_math, test_l1c_daxpy_z2);
  tcase_add_test(tc_l1c_math, test_l1c_dxmuly_z);
  tcase_add_test(tc_l1c_math, test_l1c_daxpy_z);
  tcase_add_test(tc_l1c_math, test_l1c_daxpby_z);
  tcase_add_test(tc_l1c_math, test_l1c_dsum);
  tcase_add_test(tc_l1c_math, test_l1c_dnorm1);
  tcase_add_test(tc_l1c_math, test_l1c_dlogsum);
  tcase_add_test(tc_l1c_math, test_l1c_dnrm2_err);
  tcase_add_test(tc_l1c_math, test_l1c_dnrm2_rel_err);
  tcase_add_test(tc_l1c_math, test_l1c_max_vec);
  tcase_add_test(tc_l1c_math, test_l1c_abs_vec);
  tcase_add_test(tc_l1c_math, l1c_math_max);
  tcase_add_test(tc_l1c_math, l1c_math_min);
  tcase_add_test(tc_l1c_math, l1c_math_imax);
  tcase_add_test(tc_l1c_math, l1c_math_imin);
  /*Add test cases to the suite */
  suite_add_tcase(s, tc_l1c_math);

  return s;

}
