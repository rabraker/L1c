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

#include "nesta.h"

/* Tolerances and things */
#include "test_constants.h"
#include  "l1c_math.h"


// static cJSON *test_data_json;

/* Defined in test_l1c.c*/
extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;


/* Initialize l1c_L1qcOpts struct. Even though for some tests
   we dont need all of them set, we should set them to clean up the ouput
   of valgrind when tracking down other problems.
*/

typedef struct NestaTestData {
  l1c_int m;
  l1c_int n;

  double mu;
  double sigma;
  double L;
  double tol;
  double fx_exp;

  l1c_int *pix_idx;
  double *xk;
  double *b;
  double *g;
  // double *yk;
  double *yk_exp;
  double *gradf_exp;

  //
  int n_continue;
  double beta_mu;
  double beta_tol;

  l1c_NestaProb *NP;

}NestaTestData;

static void free_generic_data(NestaTestData *Tdat);



static void init_generic_data(NestaTestData *dat){

  /* !!!!!!!!!!!!!!!! LOAD THIS FROM JSON LATER !!!!!!!!!!!!!!!!!!!*/
  dat->n_continue = 5;
  dat->tol = 1e-3;
  dat->beta_mu = 0;
  dat->beta_tol = 0;

  char *fpath_generic = fullfile(test_data_dir, "nesta_data.json");
  cJSON *json_data;


  int m=0, n=0, status=0;
  l1c_AxFuns ax_funs;

  if (load_file_to_json(fpath_generic, &json_data)){
    fprintf(stderr, "Error loading data in %s of %s\n", __func__, __FILE__);
    ck_abort();
  }
  status +=extract_json_double_array(json_data, "xk", &dat->xk, &m);
  status +=extract_json_double_array(json_data, "b", &dat->b, &n);
  status +=extract_json_double_array(json_data, "g", &dat->g, &m);
  status +=extract_json_double_array(json_data, "yk_exp", &dat->yk_exp, &m);
  status +=extract_json_double_array(json_data, "gradf_exp", &dat->gradf_exp, &m);
  status +=extract_json_int_array(json_data, "pix_idx", &dat->pix_idx, &n);
  status +=extract_json_double(json_data, "fx_exp", &dat->fx_exp);
  status +=extract_json_double(json_data, "mu", &dat->mu);
  status +=extract_json_double(json_data, "L", &dat->L);
  status +=extract_json_double(json_data, "sigma", &dat->sigma);

  if (status){
    fprintf(stderr, "Failed to load json data in %s\n", __func__);
    goto exit;
  }

  dat->n = n;
  dat->m = m;

  if(l1c_dct1_setup(dat->m, dat->n, dat->pix_idx, &ax_funs) ){
    fprintf(stderr, "Failed to initialize DCT in %s\n", __func__);
    status += 1;
    goto exit;
  }

  dat->NP = _l1c_NestaProb_new(n, m);

  if (!dat->NP){
    fprintf(stderr, "Failed to initialize Nesta Problem (in %s)\n", __func__);
    status += 1;
    goto exit;
  }


  l1c_nesta_setup(dat->NP, &dat->beta_mu, &dat->beta_tol, dat->n_continue, dat->b, ax_funs, dat->sigma,
                  dat->mu, dat->tol, dat->L, L1C_ANALYSIS);
  /* Copy the test data into NP*/
  cblas_dcopy(m, dat->xk, 1, dat->NP->xk, 1);
  cblas_dcopy(n, dat->b, 1, dat->NP->b, 1);

 exit:
  free(fpath_generic);
  cJSON_Delete(json_data);

  if (status){
    fprintf(stderr, "Error initializing data\n");
    free_generic_data(dat);
    ck_abort();
  }

}

static void free_generic_data(NestaTestData *Tdat){

    l1c_free_double(Tdat->b);
    l1c_free_double(Tdat->xk);

    l1c_free_double(Tdat->g);
    l1c_free_double(Tdat->yk_exp);
    l1c_free_double(Tdat->gradf_exp);

    free(Tdat->pix_idx);

    if (Tdat->NP){
      Tdat->NP->ax_funs.destroy();
    }

    l1c_free_nesta_problem(Tdat->NP);

}


START_TEST (test_nesta_project)
{

  double *yk=NULL;
  NestaTestData Tdat;
  init_generic_data(&Tdat);

  yk = l1c_calloc_double(Tdat.m);


  l1c_nesta_project(Tdat.NP, Tdat.xk, Tdat.g, yk);

  ck_assert_double_array_eq_tol(Tdat.m, yk, Tdat.yk_exp, TOL_DOUBLE);

  free_generic_data(&Tdat);

}END_TEST


START_TEST (test_nesta_feval)
{

  NestaTestData Tdat;
  init_generic_data(&Tdat);


  l1c_nesta_feval(Tdat.NP);

  ck_assert_double_array_eq_tol(Tdat.m, Tdat.gradf_exp, Tdat.NP->gradf, TOL_DOUBLE);

  free_generic_data(&Tdat);

}END_TEST

START_TEST (test_l1c_new_fmean_fifo)
{
  struct l1c_fmean_fifo fifo = _l1c_new_fmean_fifo();

  /*Just make sure this doesnt segfault.*/
  fifo.f_vals[L1C_NESTA_NMEAN-1] = 1;

  ck_assert_ptr_eq(fifo.f_vals, fifo.next);

  free(fifo.f_vals);
}END_TEST

START_TEST (test_l1c_push_fmeans_fifo)
{
  int i=0;
  double val;
  struct l1c_fmean_fifo fifo = _l1c_new_fmean_fifo();

  for (i=0; i<L1C_NESTA_NMEAN+1; i++){
    val = (double)i;
    _l1c_push_fmeans_fifo(&fifo, val);
  }

  ck_assert_int_eq(L1C_NESTA_NMEAN, fifo.n_total);

  ck_assert_double_eq(fifo.f_vals[0], val);

  free(fifo.f_vals);
}END_TEST


START_TEST (test_l1c_mean_fmean_fifo)
{
  int i=0;
  double val;
  double fbar_exp = 0, fbar=0;
  struct l1c_fmean_fifo fifo = _l1c_new_fmean_fifo();

  /* Ensure this works right for n_total < L1C_NESTA_NMEAN*/
  for (i=0; i<4; i++){
    val = (double)i;
    _l1c_push_fmeans_fifo(&fifo, val);
  }
  fbar = _l1c_mean_fmean_fifo(&fifo);
  /*mean([0, 1, 2, 3]) = 6/4 = 1.5*/
  fbar_exp = 1.5;
  ck_assert_double_eq_tol(fbar, fbar_exp, TOL_DOUBLE);


  /* ----------------------------------------------------
    Ensure this works right for n_total = L1C_NESTA_NMEAN
  */
  l1c_init_vec(L1C_NESTA_NMEAN, fifo.f_vals, 0);
  fifo.next = fifo.f_vals;

  for (i=0; i<L1C_NESTA_NMEAN; i++){
    val = (double)i;
    _l1c_push_fmeans_fifo(&fifo, val);
  }
  fbar = _l1c_mean_fmean_fifo(&fifo);
  /*mean([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]) = 45/10 = 4.5 */
  fbar_exp = 4.5;
  ck_assert_double_eq_tol(fbar, fbar_exp, TOL_DOUBLE);

  /* ----------------------------------------------------
     Ensure this works right for n_total == L1C_NESTA_NMEAN
     and when we have circled around.
  */
  l1c_init_vec(L1C_NESTA_NMEAN, fifo.f_vals, 0);
  fifo.next = fifo.f_vals;

  for (i=0; i<L1C_NESTA_NMEAN+2; i++){
    val = (double)i;
    _l1c_push_fmeans_fifo(&fifo, val);
  }
  fbar = _l1c_mean_fmean_fifo(&fifo);
  /*mean([10, 11, 2, 3, 4, 5, 6, 7, 8, 9]) = 65/10 = 6.5 */
  fbar_exp = 6.5;
  ck_assert_double_eq_tol(fbar, fbar_exp, TOL_DOUBLE);


  free(fifo.f_vals);
}END_TEST


Suite *l1c_nesta_suite(void)
{
  Suite *s;

  TCase *tc_nesta, *tc_fifo;

  s = suite_create("nesta");

  tc_nesta = tcase_create("nesta");
  tcase_add_test(tc_nesta, test_nesta_project);
  tcase_add_test(tc_nesta, test_nesta_feval);

  tc_fifo = tcase_create("nesta_fifo");
  tcase_add_test(tc_fifo, test_l1c_new_fmean_fifo);
  tcase_add_test(tc_fifo, test_l1c_push_fmeans_fifo);
  tcase_add_test(tc_fifo, test_l1c_mean_fmean_fifo);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_nesta);
  suite_add_tcase(s, tc_fifo);

  return s;

}
