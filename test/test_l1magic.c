/*
This is a test suite for the

This test suite uses the libcheck framework. On my computer, this got installed
into /usr/local/lib, which was not by default found by the system. Thus, I have
to do
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib


--Reference for multiple test suites:
  https://libcheck.github.io/check/doc/check_html/check_4.html#Multiple-Suites-in-one-SRunner

-- To without forking the address space:
  1. switch comments on fstat in main()
  2. in command line, run
      #: export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:/home/arnold/gradschool/sysID/sweptSines_lib"
  3. use gdb as normal
 */


#include <stdlib.h>
#include<stdio.h>
#include <check.h>

#include "TEST.h"

int main(void)
{
  enum fork_status fstat = CK_NOFORK;
  // enum fork_status fstat = CK_FORK;
  int number_failed;
  Suite *s;
  SRunner *sr;

  //  s = cgsolve_suite();
  s = l1qc_newton_suite();
  //s = dct_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, fstat);

  // Only need to do this for ADDITIONAL suites. Otherwise it will run twice.
  srunner_add_suite(sr, cgsolve_suite());
  // srunner_add_suite(sr, l1qc_newton_suite());
  srunner_add_suite(sr, dct_suite());
  /* Run the tests */
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return(number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

