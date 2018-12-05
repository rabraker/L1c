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



#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#include "cJSON.h"
#include "json_utils.h"
#include "l1qc_newton.h"


/* Tolerances and things */
#include "test_constants.h"
// #include "test_data_ss_ff.h"

static cJSON *test_data_json;


// static int load_EMx_data(int *Nx0, double **x0, int *Nx1, double **x1, int *Nidx,
//                          int **pix_idx, char *fpath){


//   if (extract_json_double_array(test_data_json, "x0", x0, Nx0) ){
//     perror("Error Loading x\n");
//     return 1;
//   }
//   if (extract_json_double_array(test_data_json, "x1", x1, Nx1) ){
//     perror("Error Loading y_exp\n");
//     return 1;
//   }

//   if (extract_json_int_array(test_data_json, "pix_idx", pix_idx, Nidx) ){
//     perror("Error Loading pix_idx \n");
//     goto end1;
//   }
//   /* Sanity check */
//   if ( (*Nx1 != *Nidx) && (*Nx0 != *Nx1)){
//     perror("Error: Array size mismatch. Aborting\n");
//     goto end2; // We allocated all, but their sizes don't match.
//   }

//   return 0;

//  end2:
//   free(*pix_idx);
//   goto end1;
//  end1:
//   free(*x1);
//   goto end0;
//  end0:
//   free(*x0);
//   return 1;

//   // end:
//   //  return 1;
// }

START_TEST(test_f_eval)
{
  double *x, *u, *r, tau, epsilon; //loaded
  double *fu1, *fu2, fe, f, *Dwork;
  double *fu1_exp, *fu2_exp, fe_exp, f_exp;

  int N, i, status=0;
  char fpath[] = "test_data/f_eval_data.json";
  load_file_to_json(fpath, &test_data_json);
  // Inputs to f_eval
  status +=extract_json_double_array(test_data_json, "x", &x, &N);
  status +=extract_json_double_array(test_data_json, "u", &u, &N);
  status +=extract_json_double_array(test_data_json, "r", &r, &N);
  status +=extract_json_double(test_data_json, "tau", &tau);
  status +=extract_json_double(test_data_json, "epsilon", &epsilon);

  // Expected outputs
  status +=extract_json_double_array(test_data_json, "fu1_exp", &fu1_exp, &N);
  status +=extract_json_double_array(test_data_json, "fu2_exp", &fu2_exp, &N);
  status +=extract_json_double(test_data_json, "f_exp", &f_exp);
  status +=extract_json_double(test_data_json, "fe_exp", &fe_exp);
  if (status){
    perror("Error Loading json data in 'test_f_eval'. Aborting\n");
    ck_abort();
  }
  Dwork = calloc(4*N, sizeof(double));
  fu1 = calloc(N, sizeof(double));
  fu2 = calloc(N, sizeof(double));

  f_eval(N, r, x, u, tau, epsilon, fu1, fu2, &fe, &f, Dwork);

  /* ----- Now check -------*/
  ck_assert_double_eq_tol(fe_exp, fe, TOL_DOUBLE);
  for (i=0; i<N; i++){
    // printf(" fu1_exp[%d]: %.11f,     fu1: %.11f\n", i, fu1_exp[i], fu1[i]);
    ck_assert_double_eq_tol(fu1_exp[i], fu1[i], TOL_DOUBLE);
    ck_assert_double_eq_tol(fu2_exp[i], fu2[i], TOL_DOUBLE);
  }

  ck_assert_double_eq_tol(f_exp, f, TOL_DOUBLE);

  /* ----------------- Cleanup --------------- */
  free(x);
  free(u);
  free(r);

  free(fu1_exp);
  free(fu2_exp);

  free(Dwork);
  free(fu1);
  free(fu2);
}
END_TEST

START_TEST(test_sum_vec)
{
  // setup_vectors_SC();
  // setup_params_SC();
  int N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double sum_exp = 21.0;
  double sum_x = sum_vec(N, x);

  ck_assert_double_eq_tol(sum_exp, sum_x, TOL_DOUBLE);
  // setup_vectors_SC();
}
END_TEST

START_TEST(test_log_vec)
{
  // setup_vectors_SC();
  // setup_params_SC();
  int i, N = 6;
  double x[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  double alpha = 3.0;
  double log_alp_x [] = {1.098612288668110,   1.791759469228055,   2.197224577336220,
                         2.484906649788000,   2.708050201102210,   2.890371757896165};
  double y[6];

  log_vec(N, alpha, x, y);

  for (i=0; i<N; i++){
    ck_assert_double_eq_tol(log_alp_x[i], y[i], TOL_DOUBLE);
  }
  // setup_vectors_SC();
}
END_TEST




Suite *l1qc_newton_suite(void)
{
  Suite *s;

  TCase *tc_core;
  s = suite_create("l1qc_newton");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_sum_vec);
  tcase_add_test(tc_core, test_log_vec);
  tcase_add_test(tc_core, test_f_eval);

  suite_add_tcase(s, tc_core);

  return s;

}
