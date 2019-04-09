#include "config.h"

#include <check.h>
#include <stdint.h>

#include "l1c.h"
#include "check_utils.h"

START_TEST(test_l1c_malloc_double){

  #define N_trial 50
  int N = 256*256;
  double *dptr_array[N_trial];
  double *dptr=NULL;
  int i;
  uintptr_t ptr_as_int;

  /* Do this a bunch of times, to try and hedge against the possibility that the allocated
   address just happens to a multiple of 64*/
  for (i=0; i<N_trial; i++){

    dptr = l1c_malloc_double(N);
    ptr_as_int = (uintptr_t) dptr;
    ck_assert_ptr_eq(0, (void*) (ptr_as_int & (DALIGN-1)) );

    dptr_array[i] = dptr;
  }

  for (i=0; i<N_trial; i++){
    l1c_free_double(dptr_array[i]);
  }



}
END_TEST

START_TEST(test_l1c_malloc_double_2D){
#define N_trial 50
  int mcol = 256*256;
  int nrow = 3;
  double **dptr_2D=NULL;
  uintptr_t ptr_as_int;

  /* Do this a bunch of times, to try and hedge against the possibility that the allocated
     address just happens to a multiple of 64*/
  for (int i=0; i<N_trial; i++){
    dptr_2D = l1c_malloc_double_2D(nrow, mcol);

    for (int k=0; k<nrow; k++){
      ptr_as_int = (uintptr_t) dptr_2D[k];
      ck_assert_ptr_eq(0, (void*) (ptr_as_int & (DALIGN-1)) );
    }

    l1c_free_double_2D(nrow, dptr_2D);
  }

}
END_TEST


START_TEST(test_l1c_calloc_double_2D){

  int N = 5;
  double x_exp[5] = {0,0,0,0,0};
  double **dptr=NULL;

  dptr = l1c_calloc_double_2D(2, N);

  for (int i=0; i<2; i++){
    ck_assert_double_array_eq(N, x_exp, dptr[i]);
  }
  l1c_free_double_2D(2, dptr);
}
END_TEST


START_TEST(test_l1c_calloc_double){

  int N = 5;
  double x_exp[5] = {0,0,0,0,0};
  double *dptr=NULL;

  dptr = l1c_calloc_double(N);


  ck_assert_double_array_eq(N, x_exp, dptr);

  l1c_free_double(dptr);
}
END_TEST


/* Add all the test cases to our suite
 */
Suite *l1c_memory_suite(void)
{
  Suite *s;

  TCase *tc_common;
  s = suite_create("l1c_memory");
  tc_common = tcase_create("l1c_memory");

  tcase_add_test(tc_common,test_l1c_malloc_double);
  tcase_add_test(tc_common, test_l1c_calloc_double);
  tcase_add_test(tc_common, test_l1c_calloc_double_2D);
  tcase_add_test(tc_common, test_l1c_malloc_double_2D);

  suite_add_tcase(s, tc_common);

  return s;

}
