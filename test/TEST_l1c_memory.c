#include "config.h"

#include <check.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "l1c_memory.h"
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

START_TEST(test_l1c_free_double_2D_when_null){
  double **dptr_2D=NULL;

  /* Nothing should happen. If it enters the loop, we should get segfault,
     and the test fails.
  */
  l1c_free_double_2D(2, dptr_2D);

  dptr_2D = l1c_malloc_double_2D(2, 10);

  /* I dont know what to do here, I guess just let the valgrind check
     complain if this didnt work right.
  */
  l1c_free_double_2D(2, dptr_2D);

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


START_TEST(test_is_aligned){

  double *x = malloc(10 * sizeof(double));
  double *y = l1c_malloc_double(10);
  bool status_x = is_aligned(x);
  bool status_y = is_aligned(y);
  ck_assert_int_eq(0, (int)status_x);
  ck_assert_int_eq(1, (int)status_y);

  free(x);
  l1c_free_double(y);
}END_TEST


START_TEST(test_daligned_offset) {
  /*
    Example: int(0x1a018f0, base=16)=27269360
    (27269360+16)/64)*64 == 27269360+16
    and sizeof(double)=8, so we should expect 2.
   */
  uintptr_t addr1 = 0x1a018f0;
  uintptr_t addr2 = 0x1a01980;

  int offset1 = next_daligned_offset((void *)addr1);
  int offset2 = next_daligned_offset((void *)addr2);

  ck_assert_int_eq(2, offset1);
  ck_assert_int_eq(0, offset2);

  for(int i=0; i<100; i++){
    double *x = malloc(i*sizeof(double));
    int offset = next_daligned_offset(x);
    bool status = is_aligned(x+offset);

    // printf("offset=%d, status=%d, x=%p\n", offset, status, (void*)x);
    ck_assert_int_eq(1, (int)status);
    ck_assert_int_lt(offset, DALIGN/sizeof(double));
    free(x);
  }

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

  tcase_add_test(tc_common, test_is_aligned);
  tcase_add_test(tc_common, test_daligned_offset);
  tcase_add_test(tc_common, test_l1c_malloc_double);
  tcase_add_test(tc_common, test_l1c_calloc_double);
  tcase_add_test(tc_common, test_l1c_calloc_double_2D);
  tcase_add_test(tc_common, test_l1c_malloc_double_2D);
  tcase_add_test(tc_common, test_l1c_free_double_2D_when_null);

  suite_add_tcase(s, tc_common);

  return s;

}
