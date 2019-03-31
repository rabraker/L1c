#include "config.h"

#include "l1c.h"

#include <check.h>
#include <stdint.h>


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


/* Add all the test cases to our suite
 */
Suite *l1c_memory_suite(void)
{
  Suite *s;

  TCase *tc_common;
  s = suite_create("l1c_memory");
  tc_common = tcase_create("l1c_memory");

  tcase_add_test(tc_common,test_l1c_malloc_double);

  suite_add_tcase(s, tc_common);

  return s;

}
