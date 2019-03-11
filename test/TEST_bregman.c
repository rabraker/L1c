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
#include "l1c_memory.h"
#include "bregman.h"
#include "check_utils.h"

extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;


cJSON *test_data_json;

typedef struct BregData {
  l1c_int NM;
  l1c_int n;
  l1c_int m;
  double *x;
  double *y;
  double *z;
  double *x_shrunk_exp;
  double *dx;
  double *dy;
  double *bx;
  double *by;
  double *f;
  double *rhs_exp;
  double *b;
  double *Hsolveb;
  double *Hessx;
  double *H_diag_exp;
  double gamma;
  double mu;
  double lam;
  char *fpath;
}BregData;

BregData *BD;

static void setup(void){
  BD = malloc(sizeof(BregData));
  BD->fpath = fullfile(test_data_dir, "bregman.json");

  cJSON *td_json;
  int setup_status=0;
  int tmp=0;
  if (load_file_to_json(BD->fpath, &td_json)){
    fprintf(stderr, "Error loading data in TV_suite\n");
    ck_abort();
  }
  free(BD->fpath);

  setup_status +=extract_json_double_array(td_json, "x", &BD->x, &BD->NM);
  setup_status +=extract_json_double_array(td_json, "y", &BD->y, &tmp);
  setup_status +=extract_json_double_array(td_json, "z", &BD->z, &tmp);
  setup_status +=extract_json_double_array(td_json, "x_shrunk", &BD->x_shrunk_exp, &tmp);
  setup_status +=extract_json_double_array(td_json, "dx", &BD->dx, &tmp);
  setup_status +=extract_json_double_array(td_json, "dy", &BD->dy, &tmp);
  setup_status +=extract_json_double_array(td_json, "bx", &BD->bx, &tmp);
  setup_status +=extract_json_double_array(td_json, "by", &BD->by, &tmp);
  setup_status +=extract_json_double_array(td_json, "f", &BD->f, &tmp);
  setup_status +=extract_json_double_array(td_json, "rhs", &BD->rhs_exp, &tmp);
  setup_status +=extract_json_double_array(td_json, "b", &BD->b, &tmp);
  setup_status +=extract_json_double_array(td_json, "Hsolveb", &BD->Hsolveb, &tmp);
  setup_status +=extract_json_double_array(td_json, "Hessx", &BD->Hessx, &tmp);
  setup_status +=extract_json_double_array(td_json, "H_diag", &BD->H_diag_exp, &tmp);

  setup_status +=extract_json_double(td_json, "gamma", &BD->gamma);
  setup_status +=extract_json_double(td_json, "mu", &BD->mu);
  setup_status +=extract_json_double(td_json, "lam", &BD->lam);
  setup_status +=extract_json_int(td_json, "m", &BD->m);
  setup_status +=extract_json_int(td_json, "n", &BD->n);

  if (setup_status){
    fprintf(stderr, "Error loading json into test data from file: \n %s. Aborting.\n",
            BD->fpath);
    ck_abort();
  }

  cJSON_Delete(td_json);
}


static void teardown(void){
  free_double(BD->x);
  free_double(BD->y);
  free_double(BD->z);
  free_double(BD->x_shrunk_exp);
  free_double(BD->dx);
  free_double(BD->dy);
  free_double(BD->bx);
  free_double(BD->by);
  free_double(BD->f);
  free_double(BD->rhs_exp);
  free_double(BD->b);
  free_double(BD->Hsolveb);
  free_double(BD->Hessx);

  free(BD);
}

START_TEST(test_breg_anis_jacobi){
  BregFuncs bfuncs = breg_get_functions();

  double *uk = malloc_double(BD->NM);
  double *dwork = malloc_double(BD->NM);
  double *D =  malloc_double(BD->NM);

  bfuncs.hess_inv_diag(BD->n, BD->m, BD->mu, BD->lam, D);

  for (int i=0; i<(BD->NM); i++){
    bfuncs.breg_anis_jacobi(BD->n, BD->m, uk, dwork, BD->b, D, BD->lam);
  }

  ck_assert_double_array_eq_tol(BD->NM, BD->Hsolveb, uk, TOL_DOUBLE);
  free_double(uk);
}
END_TEST


START_TEST(test_breg_anis_guass_seidel){
  BregFuncs bfuncs = breg_get_functions();

  double *uk = malloc_double(BD->NM);
  cblas_dcopy(BD->NM, BD->b, 1, uk, 1);

  for (int i=0; i<(BD->NM); i++){
    bfuncs.breg_anis_guass_seidel(BD->n, BD->m, uk, BD->b, BD->mu, BD->lam);
  }

  ck_assert_double_array_eq_tol(BD->NM, BD->Hsolveb, uk, TOL_DOUBLE);
  free_double(uk);
}
END_TEST

START_TEST(test_breg_shrink1){

  BregFuncs bfuncs = breg_get_functions();

  double *x_shrunk=NULL;
  l1c_int N=0;

  N = BD->NM;
  x_shrunk = malloc_double(N);
  if ( (!x_shrunk)){
    fprintf(stderr, "error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    x_shrunk[i] =0;
  }

  bfuncs.breg_shrink1(N, BD->x, x_shrunk, BD->gamma);

  ck_assert_double_array_eq_tol(N, BD->x_shrunk_exp, x_shrunk, TOL_DOUBLE_SUPER*10);

  free_double(x_shrunk);

#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST

START_TEST(test_breg_rhs){
  BregFuncs bfuncs = breg_get_functions();
  double *dwork1, *dwork2, *rhs;

  rhs = malloc_double(BD->NM);
  dwork1 = malloc_double(BD->NM);
  dwork2 = malloc_double(BD->NM);
  if ( (!rhs || !dwork1 || !dwork2) ){
    fprintf(stderr, "error allocating memory\n");
    ck_abort();
  }
  bfuncs.breg_anis_rhs(BD->n, BD->m, BD->f, BD->dx, BD->bx, BD->dy,
           BD->by, rhs, BD->mu, BD->lam, dwork1, dwork2);
  ck_assert_double_array_eq_tol(BD->NM, BD->rhs_exp, rhs, TOL_DOUBLE_SUPER*10);


  free_double(rhs);
  free_double(dwork1);
  free_double(dwork2);

#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_breg_hess_inv_diag){
  BregFuncs bfuncs = breg_get_functions();
  double *D;

  D = malloc_double(BD->NM);
  if ( (!D) ){
    fprintf(stderr, "error allocating memory\n");
    ck_abort();
  }

  bfuncs.hess_inv_diag(BD->n, BD->m, BD->mu, BD->lam, D);

  ck_assert_double_array_eq_tol(BD->NM, BD->H_diag_exp, D, TOL_DOUBLE_SUPER*10);


  free_double(D);


#ifdef _USEMKL_
  mkl_free_buffers();
#endif

}
END_TEST


START_TEST(test_breg_mxpy_z){

  BregFuncs bfuncs = breg_get_functions();
  double *z;
  l1c_int N=BD->NM;


  z = malloc_double(N);
  if ( (!z) ){
    fprintf(stderr, "error allocating memory\n");
    ck_abort();
  }
  for (int i=0; i<N; i++){
    z[i] =1;
  }
  bfuncs.breg_mxpy_z(N, BD->x, BD->y, z);

  ck_assert_double_array_eq_tol(N, BD->z, z, TOL_DOUBLE_SUPER*10);


  free_double(z);

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

  TCase *tc_breg;
  s = suite_create("bregman");
  tc_breg = tcase_create("breg");
  tcase_add_checked_fixture(tc_breg, setup, teardown);

  tcase_add_test(tc_breg, test_breg_shrink1);
  tcase_add_test(tc_breg, test_breg_mxpy_z);
  tcase_add_test(tc_breg, test_breg_hess_inv_diag);
  tcase_add_test(tc_breg, test_breg_rhs);
  tcase_add_test(tc_breg, test_breg_anis_guass_seidel);
  tcase_add_test(tc_breg, test_breg_anis_jacobi);
  suite_add_tcase(s, tc_breg);

  return s;

}
