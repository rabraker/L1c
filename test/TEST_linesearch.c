/*
This is a test suite for the l1qc_newton library.

 */
#include "config.h"

#define CK_FLOATING_DIG 17

#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <cjson/cJSON.h>

#include "cblas.h"
#include "l1c.h"
#include "check_utils.h"
#include "json_utils.h"
#include "l1c_timing.h"

#include "l1qc_newton.h"

/* Tolerances and things */
#include "test_constants.h"
#include  "l1c_math.h"
#include "linesearch.h"



// static cJSON *test_data_json;

/* Defined in test_l1c.c*/
extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;

struct quad_funx_data {
  double b; // affine term
  double g; // linear term
  double H; // Hessian term
};

struct quad_funxu_data {
  double b; // affine term
  double g1; // linear term
  double g2; // linear term
  double h11; // Hessian term
  double h12; // Hessian term
  double h22; // Hessian term
};


static double quad_funx(void *data, double *x){

  struct quad_funx_data *qd = (struct quad_funx_data*)data;
  double x1 = x[0];
  double H = qd->H;
  double g = qd->g;
  double b = qd->b;

  double f = 0.5 * (x1*x1)*H  + (x1)*g + b;
  return f;

}

static double quad_funxu(void *data, double *x_, double *u_){

  struct quad_funxu_data *qd = (struct quad_funxu_data*)data;

  double h11 = qd->h11;
  double h12 = qd->h12;
  double h22 = qd->h22;

  double g1 = qd->g1;
  double g2 = qd->g2;
  double b  = qd->b;
  double x = x_[0];
  double u = u_[0];

  double q = 0.5 * ( h11*x*x + 2*h12*x*u + h22*u*u);
  double g = g1 * x + g2 * u;

  return q + g + b;
}

void quad_funxu_dxu(void *data, double *xu, double *gradf, double *dxu){
  /* Compute the gradient and descent direction for a toy quadratic function.*/
  struct quad_funxu_data *qd = (struct quad_funxu_data*)data;

  double h11 = qd->h11;
  double h12 = qd->h12;
  double h22 = qd->h22;

  double g1 = qd->g1;
  double g2 = qd->g2;
  double x = xu[0];
  double u = xu[1];

  gradf[0] = x*h11 + u*h12 + g1;
  gradf[1] = x*h12 + u*h22 + g2;

  /* dxdu = -H^{-1}
     gradf
     H^{-1} = [h22  -h12;
              -h12   h11]
  */
  double det_H = h11*h22 - h12*h12;

  dxu[0] = - (h22 * gradf[0] - h12*gradf[1])/det_H;
  dxu[1] = - (-h12 * gradf[0] + h11 *gradf[1])/det_H;
}


START_TEST(test_linesearch_xu)
{
  /* Check the single (x) parameter linesearch. */
  LSStat ls_stat;
  LSParams ls_params;

  double *xu=NULL, *dxu=NULL, *gradf=NULL, **dwork2=NULL;
  double  fx0=0;
  double xp_exp=0, up_exp=0;
  double fp_exp;
  double flin_exp=0, s_exp=0, tmp=0;
  double g_dot_dx=0;
  double h11=3, h12=-1, h22=10, g1=1.5, g2=4, b=2;
  int n = 1, ls_iter_exp=0;

  /* Even though our test function is scalar, we need this to ensure alignment.*/
  xu = l1c_calloc_double(2);
  gradf = l1c_calloc_double(2);
  dxu = l1c_calloc_double(2);

  dwork2 = l1c_calloc_double_2D(2, 1);

  ls_params = (LSParams){.alpha = 0.6,
                         .beta = 0.9,
                         .s = 1.0};

  struct quad_funxu_data qd = {.h11=h11, .h12=h12, .h22=h22, .g1=g1, .g2=g2, .b=b};
  xu[0] = -1;
  xu[1] = 2;

  fx0 = quad_funxu((void*)(&qd), xu, xu+1);
  quad_funxu_dxu((void*)(&qd), xu, gradf, dxu);

  g_dot_dx = gradf[0] * dxu[0] + gradf[1] * dxu[1];
  printf("fx0 = %f, gradf[0] = %f, gradf[1] = %f\n", fx0, gradf[0], gradf[1]);
  printf("dx=%f, du=%f, gdotdxu=%f\n", dxu[0], dxu[1], g_dot_dx);
  /*
    Because we have a pure quadratic, we can predict how many iterations of the linesearch
    should happen. We should have s <= beta^k, where k = ceil(log(2*(1-alp))/log(beta));
  */
  tmp = log(2 * (1 - ls_params.alpha))/log(ls_params.beta);
  /* +1 because k starts from zero. */
  ls_iter_exp = (int) ceil(tmp) + 1;

  s_exp = pow(ls_params.beta, (double)(ls_iter_exp - 1));

  xp_exp = xu[0] + s_exp * dxu[0];
  up_exp = xu[1] + s_exp * dxu[1];
  fp_exp = quad_funxu((void*)(&qd), &xp_exp, &up_exp);
  flin_exp = fx0 + ls_params.alpha * s_exp * g_dot_dx;

  /* --------------- Compute ------------------------- */
  ls_stat = l1c_linesearch_xu(n, xu, xu+1, dxu, dxu+1, &fx0, g_dot_dx,
                              (void*)(&qd), quad_funxu, ls_params, dwork2);


  /* ------------------- Now check --------------------*/
  ck_assert_int_eq(ls_iter_exp, ls_stat.iter);
  ck_assert_double_eq_tol(s_exp, ls_stat.step, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(n, (&xp_exp), xu, TOL_DOUBLE);

  ck_assert_double_eq_tol(fp_exp, fx0, 1e-7);
  ck_assert_double_eq_tol(flin_exp, ls_stat.flin, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */
  l1c_free_double(xu);
  l1c_free_double(gradf);
  l1c_free_double(dxu);
  l1c_free_double_2D(2, dwork2);

}
END_TEST


START_TEST(test_linesearch)
{
  /* Check the single (x) parameter linesearch. */
  LSStat ls_stat;
  LSParams ls_params;

  double *x=NULL, *dx=NULL, gradf=0, *dwork=NULL;
  double  fx0=0;
  double *xp_exp=NULL;
  double fp_exp;
  double flin_exp=0, s_exp=0, tmp=0;
  double g_dot_dx=0;
  int n = 1, ls_iter_exp=0;

  /* Even though our test function is scalar, we need this to ensure alignment.*/
  x = l1c_calloc_double(1);
  xp_exp = l1c_calloc_double(1);
  dx = l1c_calloc_double(1);
  dwork = l1c_calloc_double(1);

  ls_params = (LSParams){.alpha = 0.6,
                         .beta = 0.9,
                         .s = 1.0};

  /* f(x) = 1/2 x*3x + 4*x + 2 */
  struct quad_funx_data qd = {.H = 3, .g = 4, .b = 2};

  x[0] = 2;
  fx0 = quad_funx((void*)(&qd), x);
  /* f'(x) = 3x + 4 */
  gradf = qd.H * x[0] + qd.g;
  /* Newton Descent direction: - H^{-1} * delF(xo) */
  dx[0] = - gradf / qd.H;
  g_dot_dx = gradf * dx[0];

  /*
    Because we have a pure quadratic, we can predict how many iterations of the linesearch
    should happen. We should have s <= beta^k, where k = ceil(log(2*(1-alp))/log(beta));
  */
  tmp = log(2 * (1 - ls_params.alpha))/log(ls_params.beta);
  /* +1 because k starts from zero. */
  ls_iter_exp = (int) ceil(tmp) + 1;

  s_exp = pow(ls_params.beta, (double)(ls_iter_exp - 1));

  xp_exp[0] = x[0] + s_exp * dx[0];
  fp_exp = quad_funx((void*)(&qd), xp_exp);
  flin_exp = fx0 + ls_params.alpha * s_exp * g_dot_dx;


  /* --------------- Compute ------------------------- */
  ls_stat = l1c_linesearch(n, x, dx, &fx0, g_dot_dx, (void*)(&qd), quad_funx, ls_params, dwork);


  /* ------------------- Now check --------------------*/
  ck_assert_int_eq(ls_iter_exp, ls_stat.iter);
  ck_assert_double_eq_tol(s_exp, ls_stat.step, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(n, xp_exp, x, TOL_DOUBLE);

  ck_assert_double_eq_tol(fp_exp, fx0, 1e-7);
  ck_assert_double_eq_tol(flin_exp, ls_stat.flin, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */
  l1c_free_double(x);
  l1c_free_double(xp_exp);
  l1c_free_double(dx);
  l1c_free_double(dwork);

}
END_TEST


START_TEST(test_linesearch_fail)
{
  /* Ensure we do the right thing when the linesearch fails. In particular, we expect
     that the linesearch returns to us the original data.
  */
  LSStat ls_stat;
  LSParams ls_params;

  double *x=NULL, *dx=NULL, gradf=0, *dwork=NULL;
  double  fx0=0;
  double g_dot_dx=0, x_orig=0;
  int n = 1;

  /* Even though our test function is scalar, we need this to ensure alignment.*/
  x = l1c_calloc_double(1);
  dx = l1c_calloc_double(1);
  dwork = l1c_calloc_double(1);

  ls_params = (LSParams){.alpha = 0.6, .beta = 0.9, .s = 1.0};

  /* f(x) = 1/2 x*3x + 4*x + 2 */
  struct quad_funx_data qd = {.H = 3, .g = 4, .b = 2};

  x[0] = 2;
  fx0 = quad_funx((void*)(&qd), x);
  /* f'(x) = 3x + 4 */
  gradf = qd.H * x[0] + qd.g;
  /* Newton Descent direction: - H^{-1} * delF(xo), but make
     dx negative, so linesearch will fail. */
  dx[0] = -(-gradf / qd.H);
  g_dot_dx = gradf * dx[0];

  x_orig = x[0];
  /*
    Because we have a pure quadratic, we can predict how many iterations of the linesearch
    should happen. We should have s <= beta^k, where k = ceil(log(2*(1-alp))/log(beta));
  */

  /* --------------- Compute ------------------------- */
  ls_stat = l1c_linesearch(n, x, dx, &fx0, g_dot_dx, (void*)(&qd), quad_funx, ls_params, dwork);

  /* ------------------- Now check --------------------*/
  ck_assert_int_eq(ls_stat.iter, MAX_LINESEARCH_ITER);
  ck_assert_double_eq_tol(pow(ls_params.beta, MAX_LINESEARCH_ITER), ls_stat.step, TOL_DOUBLE);
  /* x should not have changed.*/
  ck_assert_double_eq_tol(x_orig, x[0], TOL_DOUBLE);


  /* ----------------- Cleanup --------------- */
  l1c_free_double(x);
  l1c_free_double(dx);
  l1c_free_double(dwork);

}
END_TEST




Suite *l1c_linesearch_suite(void)
{
  Suite *s;

  TCase *tc_linesearch;

  s = suite_create("l1c_linesearch");

  tc_linesearch = tcase_create("linesearch");
  tcase_add_test(tc_linesearch, test_linesearch);
  tcase_add_test(tc_linesearch, test_linesearch_fail);
  tcase_add_test(tc_linesearch, test_linesearch_xu);


  /*Add test cases to the suite */
  suite_add_tcase(s, tc_linesearch);

  return s;

}
