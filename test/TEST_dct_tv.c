/*
  Tests for the 1D DCT-based transforms.

 */
#include "config.h"


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>
#include <cjson/cJSON.h>


/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "json_utils.h"
#include "check_utils.h"
#include "l1c.h"

#define TOL_LARGE_DCT TOL_DOUBLE_SUPER

extern char* fullfile(char *base_path, char *name);
extern char *test_data_dir;

static void check_ax_fun_properties();

typedef struct DctData{
  l1c_int *pix_idx;
  double *MtRty_RMx_exp;
  double *MtRty_RMx_act;
  double *MtRty_act;
  double *MtRty_exp;
  double *RMx_act;
  double *RMx_exp;

  double *Mx_act;
  double *Mx_exp;
  double *Mty_act;
  double *Mty_exp;

  double *Rx_act;
  double *Rx_exp;
  double *Rty_act;
  double *Rty_exp;

  double *x_in;
  double *y_in;
  double *z_in;

  double alp_v;
  double alp_h;
  /* Transform is n by mtot. mtot = mrow*mcol.*/
  l1c_int mrow;
  l1c_int mcol;
  l1c_int mtot;
  l1c_int p;
  l1c_int n;
  BpMode bp_mode;
  DctMode dct_mode;
  int setup_status;
  char *fpath;

}DctData;

/* Global variable for all test cases.
 */

static DctData *dctd;
static l1c_AxFuns ax_funs;

/* We the tcase_add_checked_fixture method. setup() and teardown are called by
   the associated setup and teardown functions for each test case. This allows us
   to use a different file path for each.
 */
// static void setup(DctData *dct_dat, char *fname){
static void setup(char *fname, BpMode bp_mode, DctMode dct_mode){
  cJSON *test_data_json;
  int setup_status=0;

  dctd = malloc(sizeof(DctData));
  if (!dctd){
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }

  dctd->fpath = fullfile(test_data_dir, fname);
  if (!dctd->fpath){
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }

  if (load_file_to_json(dctd->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_dct\n");
    ck_abort();
  }

  dctd->bp_mode = bp_mode;
  dctd->dct_mode = dct_mode;

  setup_status +=extract_json_double_array(test_data_json, "x_in", &dctd->x_in, &dctd->mtot);
  setup_status +=extract_json_double_array(test_data_json, "y_in", &dctd->y_in, &dctd->n);
  setup_status +=extract_json_double_array(test_data_json, "z_in", &dctd->z_in, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "MtEt_EMx", &dctd->MtRty_RMx_exp, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "MtEty", &dctd->MtRty_exp, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "EMx", &dctd->RMx_exp, &dctd->n);

  setup_status +=extract_json_double_array(test_data_json, "Mx", &dctd->Mx_exp, &dctd->mtot);
  setup_status +=extract_json_double_array(test_data_json, "Mty", &dctd->Mty_exp, &dctd->p);
  setup_status +=extract_json_double_array(test_data_json, "Ex", &dctd->Rx_exp, &dctd->n);
  setup_status +=extract_json_double_array(test_data_json, "Ety", &dctd->Rty_exp, &dctd->mtot);

  setup_status +=extract_json_double(test_data_json, "alp_v", &dctd->alp_v);
  setup_status +=extract_json_double(test_data_json, "alp_h", &dctd->alp_h);

  setup_status +=extract_json_int_array(test_data_json, "pix_idx", &dctd->pix_idx, &dctd->n);
  setup_status +=extract_json_int(test_data_json, "mrow", &dctd->mrow);
  setup_status +=extract_json_int(test_data_json, "mcol", &dctd->mcol);
  setup_status +=extract_json_int(test_data_json, "p", &dctd->p);

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtRty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(dctd->mrow * dctd->mcol, dctd->mtot);

  dctd->MtRty_RMx_act = l1c_malloc_double(dctd->p);
  dctd->MtRty_act = l1c_malloc_double(dctd->p);
  dctd->RMx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mty_act = l1c_malloc_double(dctd->p);
  dctd->Rx_act = l1c_calloc_double(dctd->n);
  dctd->Rty_act = l1c_malloc_double(dctd->mtot);

  int status_setup = l1c_setup_dctTV_transforms(
          dctd->n, dctd->mrow, dctd->mcol, dctd->alp_v, dctd->alp_h, dct_mode,
          bp_mode, dctd->pix_idx, &ax_funs);
  ck_assert_int_eq(status_setup, L1C_SUCCESS);

  check_ax_fun_properties();

  cJSON_Delete(test_data_json);
}

static void destroy_dctdata(){
  l1c_free_double(dctd->x_in);
  l1c_free_double(dctd->y_in);
  l1c_free_double(dctd->z_in);

  l1c_free_double(dctd->MtRty_RMx_exp);
  l1c_free_double(dctd->MtRty_exp);
  l1c_free_double(dctd->RMx_exp);
  l1c_free_double(dctd->Mx_exp);
  l1c_free_double(dctd->Mty_exp);
  l1c_free_double(dctd->Rx_exp);
  l1c_free_double(dctd->Rty_exp);

  l1c_free_double(dctd->MtRty_RMx_act);
  l1c_free_double(dctd->MtRty_act);
  l1c_free_double(dctd->RMx_act);
  l1c_free_double(dctd->Mx_act);
  l1c_free_double(dctd->Mty_act);
  l1c_free_double(dctd->Rx_act);
  l1c_free_double(dctd->Rty_act);

  free(dctd->pix_idx);
  free(dctd->fpath);

}

static void teardown(void){
  ax_funs.destroy();
  destroy_dctdata();
  free(dctd);
}


static void setup_dct2_tv_square(void){
  BpMode bp_mode = synthesis;
  DctMode dct_mode = dct2;
  setup("dct2_tv_square.json", bp_mode, dct_mode);
}


static void setup_dct2_tv_vh_square(void){
  BpMode bp_mode = synthesis;
  DctMode dct_mode = dct2;
  setup("dct2_tv_vh_square.json", bp_mode, dct_mode);
}

static void setup_dct2_tv_v_square(void){
  BpMode bp_mode = synthesis;
  DctMode dct_mode = dct2;
  setup("dct2_tv_v_square.json", bp_mode, dct_mode);
}

static void setup_dct2_tv_h_square(void){
  BpMode bp_mode = synthesis;
  DctMode dct_mode = dct2;
  setup("dct2_tv_h_square.json", bp_mode, dct_mode);
}

static void setup_analysis(void){
  BpMode bp_mode = analysis;
  DctMode dct_mode = dct2;
  setup("dct2_tv_vh_square.json", bp_mode, dct_mode);
}


static void check_ax_fun_properties() {
    ck_assert_int_eq(ax_funs.n, dctd->n);

    ck_assert_int_eq(ax_funs.q, dctd->mtot);
    if(dctd->bp_mode == synthesis){
      ck_assert_int_eq(ax_funs.m, dctd->p);
      ck_assert_int_eq(ax_funs.p, dctd->mtot);
    }else{
      ck_assert_int_eq(ax_funs.m, dctd->mtot);
      ck_assert_int_eq(ax_funs.p, dctd->p);
    }

    ck_assert_ptr_eq(ax_funs.data, NULL);
}



START_TEST(test_MtRty)
{
  ax_funs.Aty(dctd->y_in, dctd->MtRty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->MtRty_exp,
                                dctd->MtRty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_RMx)
{
  ax_funs.Ax(dctd->z_in, dctd->RMx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->RMx_exp,
                                dctd->RMx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Mz)
{
  ax_funs.Mx(dctd->z_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Mx_exp,
                                dctd->Mx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Mtx)
{
  /* --------- Wt ---------------------*/
  ax_funs.Mty(dctd->x_in, dctd->Mty_act);

  ck_assert_double_array_eq_tol(dctd->p, dctd->Mty_exp,
                                dctd->Mty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Wz_synth) {
  /* Identity in synthesis mode.*/
  ax_funs.Wz(dctd->z_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->z_in, dctd->Mx_act,
                                TOL_DOUBLE_SUPER);
}
END_TEST

START_TEST(test_Wtx_synth) {
  /* Identity in synthesis mode.*/
  ax_funs.Wtx(dctd->z_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->z_in, dctd->Mx_act,
                                TOL_DOUBLE_SUPER);
}
END_TEST


START_TEST(test_Rx)
{
  ax_funs.Rx(dctd->x_in, dctd->Rx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Rx_exp,
                                dctd->Rx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Rty)
{
  ax_funs.Rty(dctd->y_in, dctd->Rty_act);
  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Rty_exp,
                                dctd->Rty_act, TOL_DOUBLE_SUPER);

}
END_TEST

/* --------------------------- ANALYSIS MODE CHECKS -------------- */
START_TEST(test_Wz_analysis)
{
  /* In Analysis mode, this should be the identity.*/
  ax_funs.Wz(dctd->z_in, dctd->Mx_act);
  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Mx_exp,
                                dctd->Mx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Wtx_analysis)
{
  /* In Analysis mode, this should be the identity.*/
  ax_funs.Wtx(dctd->x_in, dctd->Mty_act);
  ck_assert_double_array_eq_tol(dctd->p, dctd->Mty_exp,
                                dctd->Mty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Mx_analysis)
{
  /* In Analysis mode, this should be the identity.*/
  ax_funs.Mx(dctd->x_in, dctd->Mx_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->x_in,
                                dctd->Mx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Mty_analysis)
{
  /* In Analysis mode, this should be the identity.*/
  ax_funs.Mty(dctd->x_in, dctd->Mty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->x_in,
                                dctd->Mty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Ax_analysis)
{
  /* In Analysis mode, this should be the same as Rx.*/
  ax_funs.Rx(dctd->x_in, dctd->Rx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Rx_exp,
                                dctd->Rx_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Aty_analysis)
{
  /* In Analysis mode, this should be the same as Rty.*/
  ax_funs.Rty(dctd->y_in, dctd->Rty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Rty_exp,
                                dctd->Rty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_normW)
{
  BpMode bp_mode = analysis;
  DctMode dct_mode = dct1;
  l1c_int status=0, n=5, mrow=10, mcol=10;
  double alp_v = 4.0;
  double alp_h = 5.0;
  double normW_exp=0;
  l1c_int pix_idx[5] = {0, 3, 4, 6, 7};

  normW_exp = sqrt(1 + 4*alp_v*alp_v + 4*alp_h*alp_h);
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h,
                                      dct_mode, bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_SUCCESS);
  ck_assert_double_eq_tol(ax_funs.norm_W, normW_exp, TOL_DOUBLE_SUPER);
  ax_funs.destroy();

  alp_h = 0;
  alp_v =  0;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h,
                                      dct_mode, bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_SUCCESS);
  normW_exp = sqrt(1 + 4*alp_v*alp_v + 4*alp_h*alp_h);
  ck_assert_double_eq_tol(ax_funs.norm_W, normW_exp, TOL_DOUBLE_SUPER);
  ax_funs.destroy();

  /*We should always get 1.0 in synthesis mode.*/
  bp_mode = synthesis;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, 0.0, 0.0,
                                      dct_mode, bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_SUCCESS);
  ck_assert_double_eq_tol(ax_funs.norm_W, 1.0, TOL_DOUBLE_SUPER);
  ax_funs.destroy();

}END_TEST

/* Check that the setup function returns an error code when we expect it to. */
START_TEST(test_input_errors)
{
  DctMode dct_mode = dct1;
  l1c_int status=0, n=5, mrow=10, mcol=10;
  double alp_v = 1.0;
  double alp_h = 1.0;
  l1c_int pix_idx[5] = {0, 3, 4, 6, 7};
  BpMode bp_mode = analysis;
  /* alp_v < 0 should throw an error. */
  alp_v = -1.0;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h,
                                      dct_mode, bp_mode, pix_idx, &ax_funs);

  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  /* alp_h < 0 should throw an error. */
  alp_v = 1.0;
  alp_h = -1.0;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                      bp_mode, pix_idx, &ax_funs);

  /* alp_v > 0, requires mrow>2, mcol>2   */
  alp_v = 1.0;
  alp_h = 0.0;
  mcol = 10;
  mrow = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                      bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  mrow = 10;
  mcol = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                      bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);


  /* alp_h > 0, requires mrow>2, mcol>2   */
  alp_v = 0.0;
  alp_h = 1.0;
  mcol = 10;
  mrow = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                      bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  mrow = 10;
  mcol = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                      bp_mode, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

}
END_TEST


/* Add all the test cases to our suite
 */
Suite *dctTV_suite(void)
{
  Suite *s;
  TCase *tc_errors, *tc_dct2_only, *tc_dct2_vh, *tc_dct2_v, *tc_dct2_h;
  TCase *tc_analysis;

  s = suite_create("dctTV");

  tc_errors = tcase_create("dct_tv_errors");
  tcase_add_test(tc_errors, test_input_errors);
  tcase_add_test(tc_errors, test_normW);
  suite_add_tcase(s, tc_errors);


  tc_dct2_only = tcase_create("dctTV_dct2_only");
  tcase_add_checked_fixture(tc_dct2_only, setup_dct2_tv_square, teardown);
  tcase_add_test(tc_dct2_only, test_Mtx);
  tcase_add_test(tc_dct2_only, test_Mz);
  tcase_add_test(tc_dct2_only, test_MtRty);
  tcase_add_test(tc_dct2_only, test_Rty);
  tcase_add_test(tc_dct2_only, test_Rx);
  tcase_add_test(tc_dct2_only, test_RMx);
  tcase_add_test(tc_dct2_only, test_Wtx_synth);
  tcase_add_test(tc_dct2_only, test_Wz_synth);
  suite_add_tcase(s, tc_dct2_only);

  tc_dct2_vh = tcase_create("dctTV_dct2_tv_vh");
  tcase_add_checked_fixture(tc_dct2_vh, setup_dct2_tv_vh_square, teardown);
  tcase_add_test(tc_dct2_vh, test_Mtx);
  tcase_add_test(tc_dct2_vh, test_Mz);
  tcase_add_test(tc_dct2_vh, test_MtRty);
  tcase_add_test(tc_dct2_vh, test_Rty);
  tcase_add_test(tc_dct2_vh, test_Rx);
  tcase_add_test(tc_dct2_vh, test_RMx);
  tcase_add_test(tc_dct2_vh, test_Wtx_synth);
  tcase_add_test(tc_dct2_vh, test_Wz_synth);

  suite_add_tcase(s, tc_dct2_vh);

  tc_dct2_v = tcase_create("dctTV_dct2_tv_v");
  tcase_add_checked_fixture(tc_dct2_v, setup_dct2_tv_v_square, teardown);
  tcase_add_test(tc_dct2_v, test_Mtx);
  tcase_add_test(tc_dct2_v, test_Mz);
  tcase_add_test(tc_dct2_v, test_MtRty);
  tcase_add_test(tc_dct2_v, test_Rty);
  tcase_add_test(tc_dct2_v, test_Rx);
  tcase_add_test(tc_dct2_v, test_RMx);
  tcase_add_test(tc_dct2_v, test_Wtx_synth);
  tcase_add_test(tc_dct2_v, test_Wz_synth);

  suite_add_tcase(s, tc_dct2_v);

  tc_dct2_h = tcase_create("dctTV_dct2_tv_h");
  tcase_add_checked_fixture(tc_dct2_h, setup_dct2_tv_h_square, teardown);
  tcase_add_test(tc_dct2_h, test_Mtx);
  tcase_add_test(tc_dct2_h, test_Mz);
  tcase_add_test(tc_dct2_h, test_MtRty);
  tcase_add_test(tc_dct2_h, test_Rty);
  tcase_add_test(tc_dct2_h, test_Rx);
  tcase_add_test(tc_dct2_h, test_RMx);
  tcase_add_test(tc_dct2_h, test_Wtx_synth);
  tcase_add_test(tc_dct2_h, test_Wz_synth);

  suite_add_tcase(s, tc_dct2_h);

  /* ------------- ANALYSIS -------------------*/
  tc_analysis = tcase_create("dctTV_dct2_analysis");
  tcase_add_checked_fixture(tc_analysis, setup_analysis, teardown);
  tcase_add_test(tc_analysis, test_Wz_analysis);
  tcase_add_test(tc_analysis, test_Wtx_analysis);
  tcase_add_test(tc_analysis, test_Mx_analysis);
  tcase_add_test(tc_analysis, test_Mty_analysis);
  tcase_add_test(tc_analysis, test_Ax_analysis);
  tcase_add_test(tc_analysis, test_Aty_analysis);
  tcase_add_test(tc_analysis, test_Rty);
  tcase_add_test(tc_analysis, test_Rx);


  suite_add_tcase(s, tc_analysis);
  /* The test check what happens when pix_idx is all ones, ie,
     the identitiy. That is, these are equivalent to just taking the dct/idct
     with not sampling.
  */



  return s;
}
