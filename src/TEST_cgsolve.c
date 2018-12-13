/*
  Tests for the conjugate gradient solver.

  Libcheck availible at
  https://libcheck.github.io/

 */



#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

// #include "test_data.h"
#include "cgsolve.h"

/* Tolerances and things */
#include "test_constants.h"
/* To read in test data */
#include "cJSON.h"
#include "json_utils.h"

/* Global variables which hold data contained in
   test_data_ss_ff.h
*/


int load_small_data(double **A, double **x, double **b, int *N,
                    int *max_iter, double *tol){

  cJSON *test_data_json;

  load_file_to_json("test_data/cgsolve_small01.json", &test_data_json);

  if( extract_json_int(test_data_json, "max_iter", max_iter) )
    return 1;
  if ( extract_json_double(test_data_json, "tol", tol) )
    return 1;

  int Nx, Nb, Na = 0;
  if (extract_json_double_array(test_data_json, "x", x, &Nx) ){
    perror("Error Loading x\n");
    return 1;
  }

  if (extract_json_double_array(test_data_json, "b", b, &Nb) ){
    perror("Error Loading b\n");
    goto end0;
  }

  if (extract_json_double_array(test_data_json, "A", A, &Na) ){
    perror("Error Loading A\n");
    goto end1;
  }

  *N = Nx;
  printf("Nx = %d\n", *N);
  // print_vec(*N, *x, "x");

  if ( (Nx != Nb) || (Nb*Nb != Na) ){
    perror("Error: Array size mismatch. Aborting\n");
    goto end2; // We allocated all, but their sizes don't match.
  }

  return 0;

 end2:
  free(*A);
  goto end1;
 end1:
  free(*b);
  goto end0;
 end0:
  free(*x);
  return 1;

}


START_TEST(test_cgsolve)
{
  double tol =0.0; //= 1e-6;
  int max_iter;
  CgParams cgp;
  CgResults cgr;

  double *A, *x, *x_exp, *b, *Dwork;
  int N, i = 0;
  if (load_small_data(&A, &x_exp, &b, &N, &max_iter, &tol)){
    ck_abort_msg("Errory Loading test data\n");
  }

  cgp.verbose = 0;
  cgp.tol = tol;
  cgp.max_iter = max_iter;

  x = calloc(N, sizeof(double));
  Dwork = calloc(N*4, sizeof(double));

 printf("max-iter = %d, tol=%f\n", max_iter, tol);
  cgsolve(x, b, N, Dwork, Ax, A, &cgr, cgp);

  for (i=0; i<N; i++){
    ck_assert_double_eq_tol(x_exp[i], x[i], TOL_DOUBLE);
  }

  free(A);
  free(x);
  free(x_exp);
  free(b);
  free(Dwork);

}
END_TEST




/* Add all the test cases to our suite
 */
Suite *cgsolve_suite(void)
{
  Suite *s;

  TCase *tc_core;
  s = suite_create("cgsolve");
  tc_core = tcase_create("Core");

  tcase_add_test(tc_core, test_cgsolve);

  suite_add_tcase(s, tc_core);

  return s;

}
