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
  double *MtEty_EMx_exp;
  double *MtEty_EMx_act;
  double *MtEty_act;
  double *MtEty_exp;
  double *EMx_act;
  double *EMx_exp;

  double *Mx_act;
  double *Mx_exp;
  double *Mty_act;
  double *Mty_exp;

  double *Ex_act;
  double *Ex_exp;
  double *Ety_act;
  double *Ety_exp;

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
static void setup(char *fname){
  cJSON *test_data_json;
  int setup_status=0;

  dctd = malloc(sizeof(DctData));
  if (!dctd){
    fprintf(stderr, "Memory allocation failed in %s\n", __func__);
    ck_abort();
  }

  dctd->fpath = fullfile(test_data_dir, fname);

  if (load_file_to_json(dctd->fpath, &test_data_json)){
    fprintf(stderr, "Error loading data in test_dct\n");
    ck_abort();
  }

  setup_status +=extract_json_double_array(test_data_json, "x_in", &dctd->x_in, &dctd->mtot);
  setup_status +=extract_json_double_array(test_data_json, "y_in", &dctd->y_in, &dctd->n);
  setup_status +=extract_json_double_array(test_data_json, "z_in", &dctd->z_in, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "MtEt_EMx", &dctd->MtEty_EMx_exp, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "MtEty", &dctd->MtEty_exp, &dctd->p);

  setup_status +=extract_json_double_array(test_data_json, "EMx", &dctd->EMx_exp, &dctd->n);

  setup_status +=extract_json_double_array(test_data_json, "Mx", &dctd->Mx_exp, &dctd->mtot);
  setup_status +=extract_json_double_array(test_data_json, "Mty", &dctd->Mty_exp, &dctd->p);
  setup_status +=extract_json_double_array(test_data_json, "Ex", &dctd->Ex_exp, &dctd->n);
  setup_status +=extract_json_double_array(test_data_json, "Ety", &dctd->Ety_exp, &dctd->mtot);

  setup_status +=extract_json_double(test_data_json, "alp_v", &dctd->alp_v);
  setup_status +=extract_json_double(test_data_json, "alp_h", &dctd->alp_h);

  setup_status +=extract_json_int_array(test_data_json, "pix_idx", &dctd->pix_idx, &dctd->n);
  setup_status +=extract_json_int(test_data_json, "mrow", &dctd->mrow);
  setup_status +=extract_json_int(test_data_json, "mcol", &dctd->mcol);
  setup_status +=extract_json_int(test_data_json, "p", &dctd->p);

  if (setup_status){
    fprintf(stderr, "Error Loading json into program data in 'test_MtEty_large()'. Aborting\n");
    ck_abort();
  }
  ck_assert_int_eq(dctd->mrow * dctd->mcol, dctd->mtot);

  dctd->MtEty_EMx_act = l1c_malloc_double(dctd->p);
  dctd->MtEty_act = l1c_malloc_double(dctd->p);
  dctd->EMx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mx_act = l1c_malloc_double(dctd->mtot);
  dctd->Mty_act = l1c_malloc_double(dctd->p);
  dctd->Ex_act = l1c_calloc_double(dctd->n);
  dctd->Ety_act = l1c_malloc_double(dctd->mtot);


  int status_setup = l1c_setup_dctTV_transforms(dctd->n, dctd->mrow, dctd->mcol,
                                                dctd->alp_v, dctd->alp_h,
                                                dctd->pix_idx, &ax_funs);
  ck_assert_int_eq(status_setup, L1C_SUCCESS);

  check_ax_fun_properties();

  cJSON_Delete(test_data_json);
}

static void destroy_dctdata(){
  l1c_free_double(dctd->x_in);
  l1c_free_double(dctd->y_in);
  l1c_free_double(dctd->z_in);

  l1c_free_double(dctd->MtEty_EMx_exp);
  l1c_free_double(dctd->MtEty_exp);
  l1c_free_double(dctd->EMx_exp);
  l1c_free_double(dctd->Mx_exp);
  l1c_free_double(dctd->Mty_exp);
  l1c_free_double(dctd->Ex_exp);
  l1c_free_double(dctd->Ety_exp);

  l1c_free_double(dctd->MtEty_EMx_act);
  l1c_free_double(dctd->MtEty_act);
  l1c_free_double(dctd->EMx_act);
  l1c_free_double(dctd->Mx_act);
  l1c_free_double(dctd->Mty_act);
  l1c_free_double(dctd->Ex_act);
  l1c_free_double(dctd->Ety_act);

  free(dctd->pix_idx);
  free(dctd->fpath);

}

static void teardown(void){
  ax_funs.destroy();
  destroy_dctdata();
  free(dctd);
}


static void setup_dct2_tv_square(void){
  setup("dct2_tv_square.json");
}


static void setup_dct2_tv_vh_square(void){
  setup("dct2_tv_vh_square.json");
}

static void setup_dct2_tv_v_square(void){
  setup("dct2_tv_v_square.json");
}

static void setup_dct2_tv_h_square(void){
  setup("dct2_tv_h_square.json");
}



static void check_ax_fun_properties() {
    ck_assert_int_eq(ax_funs.n, dctd->n);
    ck_assert_int_eq(ax_funs.p, dctd->p);
    ck_assert_int_eq(ax_funs.m, dctd->mrow * dctd->mcol);

    ck_assert_ptr_eq(ax_funs.data, NULL);
}



START_TEST(test_MtEty)
{
  ax_funs.Aty(dctd->y_in, dctd->MtEty_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->MtEty_exp,
                                dctd->MtEty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_EMx)
{
  ax_funs.Ax(dctd->z_in, dctd->EMx_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->EMx_exp,
                                dctd->EMx_act, TOL_DOUBLE_SUPER);

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

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Mty_exp,
                                dctd->Mty_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Ex)
{
  ax_funs.Ex(dctd->x_in, dctd->Ex_act);

  ck_assert_double_array_eq_tol(dctd->n, dctd->Ex_exp,
                                dctd->Ex_act, TOL_DOUBLE_SUPER);

}
END_TEST

START_TEST(test_Ety)
{
  ax_funs.Ety(dctd->y_in, dctd->Ety_act);

  ck_assert_double_array_eq_tol(dctd->mtot, dctd->Ety_exp,
                                dctd->Ety_act, TOL_DOUBLE_SUPER);

}
END_TEST


/* Check that the setup function returns an error code when we expect it to. */
START_TEST(test_input_errors)
{
  l1c_int status=0, n=5, mrow=10, mcol=10;
  double alp_v = 1.0;
  double alp_h = 1.0;
  l1c_int pix_idx[5] = {0, 3, 4, 6, 7};

  /* alp_v < 0 should throw an error. */
  alp_v = -1.0;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);

  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  /* alp_h < 0 should throw an error. */
  alp_v = 1.0;
  alp_h = -1.0;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);

  /* alp_v > 0, requires mrow>2, mcol>2   */
  alp_v = 1.0;
  alp_h = 0.0;
  mcol = 10;
  mrow = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  mrow = 10;
  mcol = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);


  /* alp_h > 0, requires mrow>2, mcol>2   */
  alp_v = 0.0;
  alp_h = 1.0;
  mcol = 10;
  mrow = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

  mrow = 10;
  mcol = 1;
  status = l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);
  ck_assert_int_eq(status, L1C_INVALID_ARGUMENT);

}
END_TEST


/* Add all the test cases to our suite
 */
Suite *dctTV_suite(void)
{
  Suite *s;
  TCase *tc_errors, *tc_dct2_only, *tc_dct2_vh, *tc_dct2_v, *tc_dct2_h;
  s = suite_create("dctTV");

  tc_errors = tcase_create("dct_tv_errors");
  tcase_add_test(tc_errors, test_input_errors);

  suite_add_tcase(s, tc_errors);


  tc_dct2_only = tcase_create("dctTV_dct2_only");
  tcase_add_checked_fixture(tc_dct2_only, setup_dct2_tv_square, teardown);
  tcase_add_test(tc_dct2_only, test_Mtx);
  tcase_add_test(tc_dct2_only, test_Mz);
  tcase_add_test(tc_dct2_only, test_MtEty);
  tcase_add_test(tc_dct2_only, test_Ety);
  tcase_add_test(tc_dct2_only, test_Ex);
  tcase_add_test(tc_dct2_only, test_EMx);

  suite_add_tcase(s, tc_dct2_only);

  tc_dct2_vh = tcase_create("dctTV_dct2_tv_vh");
  tcase_add_checked_fixture(tc_dct2_vh, setup_dct2_tv_vh_square, teardown);
  tcase_add_test(tc_dct2_vh, test_Mtx);
  tcase_add_test(tc_dct2_vh, test_Mz);
  tcase_add_test(tc_dct2_vh, test_MtEty);
  tcase_add_test(tc_dct2_vh, test_Ety);
  tcase_add_test(tc_dct2_vh, test_Ex);
  tcase_add_test(tc_dct2_vh, test_EMx);

  suite_add_tcase(s, tc_dct2_vh);

  tc_dct2_v = tcase_create("dctTV_dct2_tv_v");
  tcase_add_checked_fixture(tc_dct2_v, setup_dct2_tv_v_square, teardown);
  tcase_add_test(tc_dct2_v, test_Mtx);
  tcase_add_test(tc_dct2_v, test_Mz);
  tcase_add_test(tc_dct2_v, test_MtEty);
  tcase_add_test(tc_dct2_v, test_Ety);
  tcase_add_test(tc_dct2_v, test_Ex);
  tcase_add_test(tc_dct2_v, test_EMx);

  suite_add_tcase(s, tc_dct2_v);

  tc_dct2_h = tcase_create("dctTV_dct2_tv_h");
  tcase_add_checked_fixture(tc_dct2_h, setup_dct2_tv_h_square, teardown);
  tcase_add_test(tc_dct2_h, test_Mtx);
  tcase_add_test(tc_dct2_h, test_Mz);
  tcase_add_test(tc_dct2_h, test_MtEty);
  tcase_add_test(tc_dct2_h, test_Ety);
  tcase_add_test(tc_dct2_h, test_Ex);
  tcase_add_test(tc_dct2_h, test_EMx);

  suite_add_tcase(s, tc_dct2_h);

  /* The test check what happens when pix_idx is all ones, ie,
     the identitiy. That is, these are equivalent to just taking the dct/idct
     with not sampling.
  */



  return s;

}
