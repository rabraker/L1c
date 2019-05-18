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

int _l1c_l1qcProb_new( l1c_l1qcProb *Prb, l1c_int m, l1c_int n, double *b,
                       l1c_L1qcOpts params, l1c_AxFuns ax_funs);
void _l1c_l1qcProb_delete(l1c_l1qcProb *Prb);

static cJSON *test_data_json;

/* Defined in test_l1c.c*/
extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;

/* Initialize l1c_L1qcOpts struct. Even though for some tests
   we dont need all of them set, we should set them to clean up the ouput
   of valgrind when tracking down other problems.
*/
void init_l1qc_opts(l1c_L1qcOpts *params){
  params->epsilon=1e-3;
  params->tau=10;
  params->mu=10;
  params->newton_tol=1e-3;
  params->newton_max_iter=50;
  params->lbiter=0;
  params->lbtol=1e-4;
  params->l1_tol=0;
  params->verbose=2;
  params->cg_verbose=0;
  params->cg_tol=1e-8;
  params->cg_maxiter=200;
  params->warm_start_cg=0;
}

typedef struct L1qcTestData {
  l1c_int m;
  l1c_int n;
  double tau0;
  double lbtol;
  double mu;
  double epsilon;
  double cgtol;
  int cgmaxiter;
  l1c_int *pix_idx;
  double *A;
  double *b;
  double *x;
  double tau;
  int lbiter;
  double *u;
  double *r;
  double f;
  double fe;
  double *fu1;
  double *fu2;
  //for smax
  double *dx_rand1;
  double *du_rand1;
  double smax;
  // For h11p
  double *sigx_rand;
  double *z_rand;
  double *r_rand;
  double *at_rrand;
  double fe_rand;
  double *h11p_z;
  // For descent data
  double *gradf;
  double *atr;
  double *ntgx;
  double *ntgu;
  double *sig11;
  double *sig12;
  double *w1p;
  double *sigx;
  double *dx;
  double *du;

  l1c_l1qcProb Prb;

}L1qcTestData;


static void init_generic_data(L1qcTestData *dat){
  char *fpath_generic = fullfile(test_data_dir, "l1qc_data.json");
  cJSON *json_data;

  int m_tmp=0, n_tmp=0, status=0;
  l1c_AxFuns ax_funs;

  if (load_file_to_json(fpath_generic, &json_data)){
    fprintf(stderr, "Error loading data in test_l1qc_newton_1iter\n");
    ck_abort();
  }
  status +=extract_json_int(json_data, "m", &dat->m);
  status +=extract_json_int(json_data, "n", &dat->n);
  status +=extract_json_double(json_data, "tau0", &dat->tau0);
  status +=extract_json_double(json_data, "lbtol", &dat->lbtol);
  status +=extract_json_double(json_data, "mu", &dat->mu);
  status +=extract_json_double(json_data, "epsilon", &dat->epsilon);
  status +=extract_json_double(json_data, "cgtol", &dat->cgtol);
  status +=extract_json_int(json_data, "cgmaxiter", &dat->cgmaxiter);
  status +=extract_json_int_array(json_data, "pix_idx", &dat->pix_idx, &m_tmp);

  status +=extract_json_double_array(json_data, "A", &dat->A, &m_tmp);
  status +=extract_json_double_array(json_data, "b", &dat->b, &n_tmp);
  status +=extract_json_double_array(json_data, "x", &dat->x, &m_tmp);

  status +=extract_json_double(json_data, "tau", &dat->tau);
  status +=extract_json_int(json_data, "lbiter", &dat->lbiter);
  status +=extract_json_double_array(json_data, "u", &dat->u, &m_tmp);
  status +=extract_json_double_array(json_data, "r", &dat->r, &n_tmp);
  status +=extract_json_double(json_data, "f", &dat->f);
  status +=extract_json_double(json_data, "fe", &dat->fe);
  status +=extract_json_double_array(json_data, "fu1", &dat->fu1, &m_tmp);
  status +=extract_json_double_array(json_data, "fu2", &dat->fu2, &m_tmp);
  //Smax
  status +=extract_json_double_array(json_data, "dx_rand1", &dat->dx_rand1, &m_tmp);
  status +=extract_json_double_array(json_data, "du_rand1", &dat->du_rand1, &m_tmp);
  status +=extract_json_double(json_data, "smax", &dat->smax);
  //h11p
  status +=extract_json_double_array(json_data, "sigx_rand", &dat->sigx_rand, &m_tmp);
  status +=extract_json_double_array(json_data, "z_rand", &dat->z_rand, &m_tmp);
  status +=extract_json_double_array(json_data, "r_rand", &dat->r_rand, &m_tmp);
  status +=extract_json_double_array(json_data, "at_rrand", &dat->at_rrand, &m_tmp);
  status +=extract_json_double(json_data, "fe_rand", &dat->fe_rand);
  status +=extract_json_double_array(json_data, "h11p_z", &dat->h11p_z, &m_tmp);
  //For descent data

  status +=extract_json_double_array(json_data, "gradf", &dat->gradf, &m_tmp);
  status +=extract_json_double_array(json_data, "atr", &dat->atr, &m_tmp);
  status +=extract_json_double_array(json_data, "ntgx", &dat->ntgx, &m_tmp);
  status +=extract_json_double_array(json_data, "ntgu", &dat->ntgu, &m_tmp);
  status +=extract_json_double_array(json_data, "sig11", &dat->sig11, &m_tmp);
  status +=extract_json_double_array(json_data, "sig12", &dat->sig12, &m_tmp);
  status +=extract_json_double_array(json_data, "w1p", &dat->w1p, &m_tmp);
  status +=extract_json_double_array(json_data, "sigx", &dat->sigx, &m_tmp);
  status +=extract_json_double_array(json_data, "dx", &dat->dx, &m_tmp);
  status +=extract_json_double_array(json_data, "du", &dat->du, &m_tmp);

  l1c_L1qcOpts params;
  init_l1qc_opts(&params);
  params.epsilon = dat->epsilon;
  params.tau = dat->tau;

  if(l1c_dct1_setup(dat->m, dat->n, dat->pix_idx, &ax_funs) ){
    ck_abort();
  }

  status += _l1c_l1qcProb_new(&(dat->Prb), dat->m, dat->n, dat->b, params, ax_funs);

  if (status){
    fprintf(stderr, "Error initializing data\n");
    ck_abort();
  }

  free(fpath_generic);
  cJSON_Delete(json_data);
}

void free_generic_data(L1qcTestData Tdat){
    free(Tdat.pix_idx);
    l1c_free_double(Tdat.A);
    l1c_free_double(Tdat.b);
    l1c_free_double(Tdat.x);

    l1c_free_double(Tdat.u);
    l1c_free_double(Tdat.r);
    l1c_free_double(Tdat.fu1);
    l1c_free_double(Tdat.fu2);
    //for smax
    l1c_free_double(Tdat.dx_rand1);
    l1c_free_double(Tdat.du_rand1);
    // For h11p
    l1c_free_double(Tdat.sigx_rand);
    l1c_free_double(Tdat.z_rand);
    l1c_free_double(Tdat.r_rand);
    l1c_free_double(Tdat.at_rrand);
    l1c_free_double(Tdat.h11p_z);
    // For descent data
    l1c_free_double(Tdat.gradf);
    l1c_free_double(Tdat.atr);
    l1c_free_double(Tdat.ntgx);
    l1c_free_double(Tdat.ntgu);
    l1c_free_double(Tdat.sig11);
    l1c_free_double(Tdat.sig12);
    l1c_free_double(Tdat.w1p);
    l1c_free_double(Tdat.sigx);
    l1c_free_double(Tdat.dx);
    l1c_free_double(Tdat.du);

    _l1c_l1qcProb_delete(&(Tdat.Prb));

    Tdat.Prb.ax_funs.destroy();
}


START_TEST (test_l1qc_check_feasible_start)
{
  l1c_L1qcOpts params;
  char *fpath_1iter = fullfile(test_data_dir, "lb_test_data_AX.json");
  double *x0=NULL, *b=NULL;
  double *A=NULL;
  int is_feasible=0;
  l1c_int m=0, n=0, mn=0,  status=0;
  l1c_AxFuns ax_funs;
  l1c_l1qcProb Prb;

  init_l1qc_opts(&params);

  if (load_file_to_json(fpath_1iter, &test_data_json)){
    fprintf(stderr, "Error loading data in test_l1qc_newton_1iter\n");
    free(fpath_1iter);
    ck_abort();
  }
  free(fpath_1iter);

  status +=extract_json_double_array(test_data_json, "x0", &x0, &m);
  status +=extract_json_double_array(test_data_json, "b", &b, &n);
  status +=extract_json_double_array(test_data_json, "A", &A, &mn);

  ck_assert_int_eq(mn, m*n);

  if (status || !A || !x0 || !b){
    status += 1;
    goto exit1;
  }
  if(l1c_setup_matrix_transforms(n, m, A, &ax_funs)){
    goto exit1;
  }

  if (_l1c_l1qcProb_new(&Prb, m, n, b, params, ax_funs) ){
    goto exit1;
  }

  Prb.epsilon = 0.0;
  is_feasible = _l1c_l1qc_check_feasible_start(&Prb, x0);
  ck_assert_int_eq(is_feasible, L1C_INFEASIBLE_START);

  Prb.epsilon = 1.0;
  is_feasible = _l1c_l1qc_check_feasible_start(&Prb, x0);
  ck_assert_int_eq(is_feasible, 0);

 exit1:
  l1c_free_double(A);
  l1c_free_double(x0);
  l1c_free_double(b);
  _l1c_l1qcProb_delete(&Prb);
  ax_funs.destroy();


}END_TEST


START_TEST (test_l1qc_newton)
{
  char *fpath_1iter = fullfile(test_data_dir, "lb_test_data_AX.json");

  double *x0=NULL, *b=NULL;
  double *x_exp=NULL, *A=NULL;
  double enrm1;

  l1c_int m=0, n=0, mn=0,  status=0;
  l1c_int T=0, TC=0;
  l1c_int *T_idx=NULL, *TC_idx=NULL;

  l1c_AxFuns ax_funs;
  l1c_LBResult lb_res;
  l1c_L1qcOpts params;

  init_l1qc_opts(&params);

  if (load_file_to_json(fpath_1iter, &test_data_json)){
    fprintf(stderr, "Error loading data in test_l1qc_newton_1iter\n");
    free(fpath_1iter);
    ck_abort();
  }
  free(fpath_1iter);

  status +=extract_json_double_array(test_data_json, "x0", &x0, &n);
  status +=extract_json_double_array(test_data_json, "x_act", &x_exp, &n);
  status +=extract_json_double_array(test_data_json, "b", &b, &m);

  status +=extract_json_double(test_data_json, "epsilon", &params.epsilon);
  status +=extract_json_double(test_data_json, "mu", &params.mu);
  status +=extract_json_double(test_data_json, "lbtol", &params.lbtol);
  status +=extract_json_double(test_data_json, "newtontol", &params.newton_tol);
  status +=extract_json_int(test_data_json, "newtonmaxiter", &params.newton_max_iter);
  status +=extract_json_double(test_data_json, "cgtol", &params.cg_tol);
  status +=extract_json_double(test_data_json, "enrm1", &enrm1);
  status +=extract_json_int(test_data_json, "cgmaxiter", &params.cg_maxiter);

  status +=extract_json_int_array(test_data_json, "T_idx", &T_idx, &T);
  status +=extract_json_int_array(test_data_json, "TC_idx", &TC_idx, &TC);
  status +=extract_json_double_array(test_data_json, "A", &A, &mn);

  ck_assert_int_eq(mn, m*n);

  double *tmp1 = l1c_calloc_double(n);
  double *tmp2 = l1c_calloc_double(n);

  if (status || !TC_idx || !T_idx || !A || !x0 || !x_exp || !b || !tmp1 || !tmp2){
    status += 1;
    goto exit1;
  }

  if(l1c_setup_matrix_transforms(m, n, A, &ax_funs)){
    goto exit1;
  }

  params.verbose = 2;
  params.warm_start_cg = 0;
  params.l1_tol = 1e-4; //wont be used.
  struct timeval tv_start, tv_end;
  tv_start = l1c_get_time();

  lb_res = l1c_l1qc_newton(n, x0, m, b, params, ax_funs);

  tv_end = l1c_get_time();
  double time_total = l1c_get_time_diff(tv_start, tv_end);
  printf("total c time: %f\n", time_total);

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
  double dnrm1_x1exp= l1c_dnorm1(n, x_exp);
  double dnrm1_xp = l1c_dnorm1(n, x0);
  ck_assert_double_le(dnrm1_xp, dnrm1_x1exp);

  /* -----------------------
    b. Check property (b) */
  ax_funs.Ax(x_exp, tmp1);
  ax_funs.Ax(x0, tmp2);
  double dnrm2_yerr = 0, yerr_k = 0;
  for(int i=0; i<n; i++){
    yerr_k = tmp1[i] - tmp2[i];
    dnrm2_yerr += yerr_k * yerr_k;
  }
  dnrm2_yerr = sqrt(dnrm2_yerr);

  ck_assert_double_le(dnrm2_yerr, params.epsilon*2);

  /* c. ------------------
     Check property (c) */
  l1c_daxpy_z(n, 1, x0, x_exp, tmp1); //h = x0 - x
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
  ck_assert_int_eq(0, lb_res.status);


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
  cJSON_Delete(test_data_json);


  if (status){
    fprintf(stderr, "Error Loading json data in 'test_l1qc_newton_1ter()'. Aborting\n");
    ck_abort();
  }

}
END_TEST


START_TEST (test_newton_init_regres1)
{

  l1c_L1qcOpts params;
  init_l1qc_opts(&params);
  l1c_int m=4;

  double x[] = {1.0, 2.0, 3.0, 4.0};
  double u[] = {0,0,0,0};
  double u_exp[] = {1.3500, 2.3000, 3.2500, 4.2000};

  params.lbtol = 0.1;
  params.mu = 0.1;

 _l1c_l1qc_newton_init(m, x, u, &params);

  ck_assert_double_array_eq_tol(m, u_exp, u,  TOL_DOUBLE_SUPER);

}
END_TEST


START_TEST (test_newton_init)
{
  L1qcTestData Tdat;
  double *u=NULL;
  int ret=0, status=0;
  l1c_L1qcOpts params;
  init_generic_data(&Tdat);
  u = l1c_malloc_double(Tdat.m);
  if ( !u ){
    fprintf(stderr, "Error Allocating Memory in 'test_newton_init()'. Aborting\n");
    status = L1C_OUT_OF_MEMORY;
    goto exit1;
  }

  params =(l1c_L1qcOpts) {.mu=Tdat.mu, .lbtol=Tdat.lbtol,
                          .epsilon=Tdat.epsilon, .lbiter = 0};

  ret= _l1c_l1qc_newton_init(Tdat.m, Tdat.x, u, &params);

  ck_assert_int_eq(0, ret);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.u, u,  TOL_DOUBLE_SUPER);
  ck_assert_double_eq_tol(Tdat.tau, params.tau,  TOL_DOUBLE_SUPER*100);
  ck_assert_int_eq(Tdat.lbiter, params.lbiter);

  params.lbiter = 1;
  ret= _l1c_l1qc_newton_init(Tdat.m, Tdat.x, u, &params);
  ck_assert_int_eq(1, params.lbiter);


 exit1:
  l1c_free_double(u);

  free_generic_data(Tdat);

  if (status){
    ck_abort();
  }

}
END_TEST


START_TEST(test_l1qc_descent_dir)
{
  int status=0;
  l1c_int m = 0;
  L1qcTestData Tdat;
  init_generic_data(&Tdat);

  l1c_l1qcProb l1qc_prob = Tdat.Prb;
  l1c_CgParams cgp = {.verbose=0, .max_iter=Tdat.cgmaxiter, .tol = Tdat.cgtol};
  l1c_CgResults cgr = {.cgres=0.0, .cgiter=0};

  m = Tdat.m;

  cblas_dcopy(Tdat.n,  Tdat.r, 1, l1qc_prob.r, 1);
  cblas_dcopy(m, Tdat.fu1, 1, l1qc_prob.fu1, 1);
  cblas_dcopy(m, Tdat.fu2, 1, l1qc_prob.fu2, 1);

  l1qc_prob.fe_val = Tdat.fe;
  l1qc_prob.tau = Tdat.tau0;

  /* We must initialize l1qc_prob.dx, because we have enabled warm starting*/
  l1c_init_vec(Tdat.m, l1qc_prob.dx, 0.0);

  _l1c_l1qc_descent_dir(&l1qc_prob, cgp, &cgr);

  /*The next three should already be checked by test_get_gradient, but we can
   do it here too, to make sure things are staying sane.*/
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.sig11, l1qc_prob.sig11, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.sig12, l1qc_prob.sig12, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.w1p, l1qc_prob.w1p, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.dx, l1qc_prob.dx, TOL_DOUBLE*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.du, l1qc_prob.du, TOL_DOUBLE*100);

  /* ----------------- Cleanup --------------- */

  free_generic_data(Tdat);

  if (status)
    ck_abort();

}
END_TEST


START_TEST(test_H11pfun)
{
  L1qcTestData Tdat;
  init_generic_data(&Tdat);

  Hess_data h11p_data;
  double *h11p_z=NULL, *z_orig=NULL;
  int status=0;

  h11p_z = l1c_malloc_double(Tdat.m);
  z_orig = l1c_malloc_double(Tdat.m);
  h11p_data.Dwork_1m = l1c_malloc_double(Tdat.m);
  if (!h11p_z || !z_orig || !h11p_data.Dwork_1m){
    fprintf(stderr, "Unable to allocate memory\n");
    status = L1C_OUT_OF_MEMORY;
    goto exit;
  }

  h11p_data.one_by_fe = 1.0/Tdat.fe_rand;
  h11p_data.one_by_fe_sqrd = 1.0/(Tdat.fe_rand * Tdat.fe_rand);
  h11p_data.atr = Tdat.at_rrand;
  h11p_data.sigx = Tdat.sigx_rand;

  h11p_data.AtAx = Tdat.Prb.ax_funs.AtAx;

  cblas_dcopy(Tdat.m, Tdat.z_rand, 1, z_orig, 1);


  _l1c_l1qc_H11pfun(Tdat.m, Tdat.z_rand, h11p_z, &h11p_data);

  ck_assert_double_array_eq_tol(Tdat.m, Tdat.h11p_z, h11p_z, TOL_DOUBLE);
  // Ensure we didnt overwrite data in z.
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.z_rand, z_orig, TOL_DOUBLE_SUPER);


  /* ----------------- Cleanup --------------- */
 exit:
  l1c_free_double(h11p_z);
  l1c_free_double(z_orig);
  l1c_free_double(h11p_data.Dwork_1m);

  free_generic_data(Tdat);

  if (status)
    ck_abort();

}
END_TEST


START_TEST(test_l1qc_hess_grad)
{
  l1c_int m = 0;
  int status = 0;

  L1qcTestData Tdat;
  init_generic_data(&Tdat);
  m = Tdat.m;

  l1c_l1qcProb l1qc_prob = Tdat.Prb;
  double *sigx = l1c_calloc_double(m);
  if (!sigx){
    status = 0;
    goto exit;
  }

  cblas_dcopy(m, Tdat.fu1, 1, l1qc_prob.fu1, 1);
  cblas_dcopy(m, Tdat.fu2, 1, l1qc_prob.fu2, 1);

  l1qc_prob.fe_val = Tdat.fe;
  l1qc_prob.tau = Tdat.tau0;

  /*-------------- Compute --------------------------- */
  _l1c_l1qc_hess_grad(&l1qc_prob, sigx, Tdat.atr);

  /* ----- Check -------*/
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.sig11, l1qc_prob.sig11, TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.sig12, l1qc_prob.sig12, TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.w1p, l1qc_prob.w1p, TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.sigx, sigx, TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.ntgu, l1qc_prob.ntgu, TOL_DOUBLE_SUPER*100);
  ck_assert_double_array_eq_tol(Tdat.m*2, Tdat.gradf, l1qc_prob.gradf, TOL_DOUBLE_SUPER*100);

  /* ----------------- Cleanup --------------- */

 exit:
  l1c_free_double(sigx);

  free_generic_data(Tdat);

  if (status)
    ck_abort();
}
END_TEST


START_TEST (test_find_max_step)
{
  L1qcTestData Tdat;
  init_generic_data(&Tdat);

  l1c_l1qcProb l1qc_prob = Tdat.Prb;
  double smax= 0.0;
  l1c_int m = Tdat.m;

  cblas_dcopy(m, Tdat.dx_rand1, 1, l1qc_prob.dx, 1);
  cblas_dcopy(m, Tdat.du_rand1, 1, l1qc_prob.du, 1);
  cblas_dcopy(m, Tdat.fu1, 1, l1qc_prob.fu1, 1);
  cblas_dcopy(m, Tdat.fu2, 1, l1qc_prob.fu2, 1);
  cblas_dcopy(Tdat.n, Tdat.r, 1, l1qc_prob.r, 1);

  l1qc_prob.epsilon = Tdat.epsilon;
  l1qc_prob.m = Tdat.m;
  l1qc_prob.n = Tdat.n;

  /* ------------- Compute ------------*/
  smax = _l1c_l1qc_find_max_step(&l1qc_prob);

  /* ------------- Check --------------*/
  ck_assert_double_eq_tol(Tdat.smax, smax,  TOL_DOUBLE);

  free_generic_data(Tdat);

}
END_TEST


START_TEST(test_f_eval)
{
  L1qcTestData Tdat;
  init_generic_data(&Tdat);

  int status = 0;
  l1c_l1qcProb l1qc_prob = Tdat.Prb;

  if (status){
    fprintf(stderr, "Error allocating memory in test_f_eval\n");
    status = L1C_OUT_OF_MEMORY;
    goto exit;
  }

  l1qc_prob.tau = Tdat.tau0;
  l1qc_prob.epsilon = Tdat.epsilon;

  /* ------------- Compute----------------------*/
  _l1c_l1qc_f_eval( &l1qc_prob, Tdat.x, Tdat.u);

  /* ------------- check -----------------------*/
  ck_assert_double_eq_tol(Tdat.fe, l1qc_prob.fe_val, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.fu1, l1qc_prob.fu1, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(Tdat.m, Tdat.fu2, l1qc_prob.fu2, TOL_DOUBLE);
  ck_assert_double_array_eq_tol(Tdat.n, Tdat.r, l1qc_prob.r, TOL_DOUBLE);

  ck_assert_double_eq_tol(Tdat.f, l1qc_prob.f_val, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */

 exit:

  free_generic_data(Tdat);

  if (status)
    ck_abort();

}
END_TEST


Suite *l1qc_newton_suite(void)
{
  Suite *s;

  TCase *tc_linesearch, *tc_newton_init, *tc_feval, *tc_gradient;
  TCase *tc_l1qc_newton, *tc_h11p;

  s = suite_create("l1qc_newton");

  tc_l1qc_newton = tcase_create("l1qc_newton");
  tcase_add_test(tc_l1qc_newton, test_l1qc_newton);

  tc_newton_init = tcase_create("newton_init");
  tcase_add_test(tc_newton_init, test_newton_init);
  tcase_add_test(tc_newton_init, test_newton_init_regres1);
  tcase_add_test(tc_newton_init, test_l1qc_check_feasible_start);

  tc_feval = tcase_create("f_eval");
  tcase_add_test(tc_feval, test_f_eval);

  tc_linesearch = tcase_create("l1qc_linesearch");
  tcase_add_test(tc_linesearch, test_find_max_step);
  // tcase_add_test(tc_linesearch, test_line_search);

  tc_gradient = tcase_create("l1qc_gradient");
  tcase_add_test(tc_gradient, test_l1qc_hess_grad);
  tcase_add_test(tc_gradient, test_l1qc_descent_dir);

  tc_h11p = tcase_create("h11p");
  tcase_add_test(tc_h11p, test_H11pfun);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_l1qc_newton);
  suite_add_tcase(s, tc_newton_init);
  suite_add_tcase(s, tc_feval);
  suite_add_tcase(s, tc_linesearch);
  suite_add_tcase(s, tc_gradient);
  suite_add_tcase(s, tc_h11p);

  return s;

}
