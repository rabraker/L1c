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
  BpMode bp_mode = analysis;
  DctMode dct_mode = dct1;
  if(l1c_setup_dctTV_transforms(dat->n, dat->m, 1, 0.0, 0.0, dct_mode,
                                bp_mode, dat->pix_idx, &ax_funs) ){
    fprintf(stderr, "Failed to initialize DCT in %s\n", __func__);
    status += 1;
    goto exit;
  }
  ax_funs.Ax = ax_funs.Rx;
  ax_funs.Aty = ax_funs.Rty;

  dat->NP = _l1c_NestaProb_new(ax_funs);

  if (!dat->NP){
    fprintf(stderr, "Failed to initialize Nesta Problem (in %s)\n", __func__);
    status += 1;
    goto exit;
  }

  l1c_NestaOpts opts = {.n_continue = dat->n_continue,
                        .sigma = dat->sigma,
                        .mu = dat->mu,
                        .tol = dat->tol,
                        .bp_mode = analysis};

  if (l1c_nesta_setup(dat->NP, &dat->beta_mu, &dat->beta_tol, dat->b, ax_funs,
                        &opts)) {
    status++;
    goto exit;
  }

  /* Copy the test data into NP*/
  cblas_dcopy(m, dat->xk, 1, dat->NP->xk, 1);
  cblas_dcopy(n, dat->b, 1, dat->NP->b, 1);

 exit:
  free(fpath_generic);
  cJSON_Delete(json_data);

  if (status){
    fprintf(stderr, "Error initializing data in %s\n", __func__);
    free_generic_data(dat);
    ck_abort();
  }

}


static void free_generic_data(NestaTestData *Tdat){

  l1c_free_double(Tdat->xk);
  l1c_free_double(Tdat->b);
  l1c_free_double(Tdat->g);
  l1c_free_double(Tdat->yk_exp);
  l1c_free_double(Tdat->gradf_exp);

  free(Tdat->pix_idx);

  if (Tdat->NP){
    Tdat->NP->ax_funs.destroy();
  }

  l1c_free_nesta_problem(Tdat->NP);

}



START_TEST (test_l1c_nesta)
{
  char *fpath_1iter = fullfile(test_data_dir, "lb_test_data_AX.json");
  cJSON *json_data;

  double *x0=NULL, *b=NULL;
  double *x_exp=NULL, *A=NULL;
  double enrm1, epsilon;

  l1c_int m=0, n=0, mn=0,  status=0;
  l1c_int T=0, TC=0;
  l1c_int *T_idx=NULL, *TC_idx=NULL;

  l1c_AxFuns ax_funs;


  if (load_file_to_json(fpath_1iter, &json_data)){
    fprintf(stderr, "Error loading data in test_l1qc_newton_1iter\n");
    free(fpath_1iter);
    ck_abort();
  }
  free(fpath_1iter);

  status +=extract_json_double_array(json_data, "x0", &x0, &m);
  status +=extract_json_double_array(json_data, "x_act", &x_exp, &m);
  status +=extract_json_double_array(json_data, "b", &b, &n);

  status +=extract_json_double(json_data, "epsilon", &epsilon);
  status +=extract_json_double(json_data, "enrm1", &enrm1);

  status +=extract_json_int_array(json_data, "T_idx", &T_idx, &T);
  status +=extract_json_int_array(json_data, "TC_idx", &TC_idx, &TC);
  status +=extract_json_double_array(json_data, "A", &A, &mn);

  ck_assert_int_eq(mn, n*m);

  double *tmp1 = l1c_calloc_double(m);
  double *tmp2 = l1c_calloc_double(m);

  if (status || !TC_idx || !T_idx || !A || !x0 || !x_exp || !b || !tmp1 || !tmp2){
    status += 1;
    goto exit1;
  }

  if(l1c_setup_matrix_transforms(n, m, A, &ax_funs)){
    goto exit0;
  }

  // matrix_transforms doesnt set normW yet. Would be better to compute in python.
  ax_funs.norm_W = 1;
  l1c_NestaOpts opts = {.mu=1e-5,
                        .sigma=epsilon,
                        .tol=1e-3,
                        .n_continue=5,
                        .bp_mode=synthesis};
  /* ------------------------------------------------------- */
  int nesta_status = l1c_nesta(m, x0, n, b, ax_funs, opts);

  ck_assert_int_eq(nesta_status, 0);
  /*
    From Stable Signal Recovery from Incomplete and Inaccurate Measurements,
    Candes, Romberg, Tao 2005, it looks like we should be able to verify
    (cir. (9)-(10)):
       a. ||x_opt||_1 <= ||x_act||_1
       b. ||Ax_opt - Ax_act||_1 <= 2*\epsilon
       c. Let h = x_opt - x_act. Let h_T = h[idx_supp], and
       h_TC = h[(1:m)!=idx_supp]. Then
       ||h_TC||_1 <= ||h_T||_1

       d. From their numerical experiments, they suggest that
          ||x_opt - x_act|| < C*eps, with C<2 (see also eq. (5)).
          with eps^2 = sigma^2(n + lambda*sqrt(2*n)). Though table 1
          shows this, I dont get the same result repeating that
          experiment, even with their software. It seems that C~=3.5
   */


  /* a. ---------------------
     Check property (a).*/
  double dnrm1_x1exp= l1c_dnorm1(m, x_exp);
  double dnrm1_xp = l1c_dnorm1(m, x0);
  ck_assert_double_le(dnrm1_xp, dnrm1_x1exp);

  /* -----------------------
    b. Check property (b) */
  ax_funs.Ax(x_exp, tmp1);
  ax_funs.Ax(x0, tmp2);
  double dnrm2_yerr = 0, yerr_k = 0;
  for(int i=0; i<m; i++){
    yerr_k = tmp1[i] - tmp2[i];
    dnrm2_yerr += yerr_k * yerr_k;
  }
  dnrm2_yerr = sqrt(dnrm2_yerr);

  ck_assert_double_le(dnrm2_yerr, epsilon*2);

  /* c. ------------------
     Check property (c) */
  l1c_daxpy_z(m, -1, x0, x_exp, tmp1); //h = x0 - x
  double h_T_nrm1=0, h_TC_nrm1=0;
  for(int i=0; i<T; i++){
    h_T_nrm1 += fabs(tmp1[T_idx[i]]);
  }
  for(int i=0; i<TC; i++){
    h_TC_nrm1 += fabs(tmp1[TC_idx[i]]);
  }
  ck_assert_double_le(h_TC_nrm1, h_T_nrm1);

  /*d. TODO: implement this check. */

  /* e. We should exit cleanly for this problem.*/
  ck_assert_int_eq(0, nesta_status);

 exit0:
  ax_funs.destroy();
 exit1:
  l1c_free_double(x0);
  l1c_free_double(x_exp);
  l1c_free_double(b);
  l1c_free_double(A);
  l1c_free_double(tmp1);
  l1c_free_double(tmp2);
  free(T_idx);
  free(TC_idx);
  cJSON_Delete(json_data);


  if (status){
    fprintf(stderr, "Error Loading json data in 'test_l1qc_newton_1ter()'. Aborting\n");
    ck_abort();
  }

}
END_TEST


START_TEST (test_nesta_project)
{

  double *yk=NULL;
  NestaTestData Tdat;
  init_generic_data(&Tdat);

  yk = l1c_calloc_double(Tdat.m);

  Tdat.NP->mu_j = Tdat.mu;
  l1c_nesta_project(Tdat.NP, Tdat.xk, Tdat.g, yk);

  ck_assert_double_array_eq_tol(Tdat.m, yk, Tdat.yk_exp, TOL_DOUBLE);

  free_generic_data(&Tdat);
  l1c_free_double(yk);

}END_TEST


START_TEST (test_nesta_feval)
{

  NestaTestData Tdat;
  init_generic_data(&Tdat);
  Tdat.NP->mu_j = Tdat.NP->mu;

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


START_TEST (test_l1c_nesta_setup)
{
  l1c_int n = 5;
  l1c_int m = 10;

  double beta_mu=0, beta_tol=0;
  double sigma = 1e-3, mu = 1e-5;
  double tol = 1e-3;

  int n_continue = 5;
  double *b = l1c_calloc_double(n);
  double *A = l1c_calloc_double(n*m);
  DctMode dct_mode = dct1;

  for (int i=0; i<n; i++){
    b[i] = ((double)rand())/(double)RAND_MAX;
  }
  for (int i=0; i<n*m; i++){
      A[i] = ((double)rand())/(double)RAND_MAX;
  }
  l1c_AxFuns ax_funs;
  l1c_setup_matrix_transforms(n, m, A, &ax_funs);

  l1c_NestaProb *NP = _l1c_NestaProb_new(ax_funs);

  l1c_NestaOpts opts = {.n_continue=5, .sigma=sigma, .mu=mu, .tol=tol, .bp_mode=analysis};
  // We are checking that setup will fail for ax_funs without analysis opertator.
  ax_funs.Wz = NULL;
  int status = l1c_nesta_setup(NP, &beta_mu, &beta_tol, b, ax_funs, &opts);

  ck_assert_int_eq(status, L1C_INCONSISTENT_ARGUMENTS);

  ax_funs.destroy();

  int pix_idx[5] = {1, 3, 4, 6, 8};
  if (l1c_setup_dct_transforms(n, m, 1, dct_mode, pix_idx, &ax_funs)){
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }


  status = l1c_nesta_setup(NP, &beta_mu, &beta_tol, b, ax_funs, &opts);

  ck_assert_int_eq(status, L1C_SUCCESS);
  ck_assert_ptr_eq(b, NP->b);
  ck_assert_int_eq(n_continue, NP->n_continue);
  ck_assert_double_eq(NP->sigma, sigma);
  ck_assert_double_eq(NP->mu, mu);
  ck_assert_double_eq(NP->tol, tol);

  double mu_j=NP->mu_j, tol_j=NP->tol_j;

  for (int i=0; i < n_continue; i ++){
    mu_j *=  beta_mu;
    tol_j *=  beta_tol;
  }

  ck_assert_double_eq_tol(mu_j, mu, TOL_DOUBLE);
  ck_assert_double_eq_tol(tol_j, tol, TOL_DOUBLE);

  l1c_free_nesta_problem(NP);
  ax_funs.destroy();
  l1c_free_double(b);
  l1c_free_double(A);

} END_TEST

Suite *l1c_nesta_suite(void)
{
  Suite *s;

  TCase *tc_nesta, *tc_fifo;

  s = suite_create("nesta");

  tc_nesta = tcase_create("nesta");
  tcase_add_test(tc_nesta, test_nesta_project);
  tcase_add_test(tc_nesta, test_nesta_feval);
  tcase_add_test(tc_nesta, test_l1c_nesta);

  tcase_add_test(tc_nesta, test_l1c_nesta_setup);

  tc_fifo = tcase_create("nesta_fifo");
  tcase_add_test(tc_fifo, test_l1c_new_fmean_fifo);
  tcase_add_test(tc_fifo, test_l1c_push_fmeans_fifo);
  tcase_add_test(tc_fifo, test_l1c_mean_fmean_fifo);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_nesta);
  suite_add_tcase(s, tc_fifo);

  return s;

}
