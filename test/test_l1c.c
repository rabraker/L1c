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

#include "config.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <check.h>

#include "TEST.h"

char *test_data_dir;
static char test_data_dir_name[] = "test_data";

char* fullfile(char *base_path, char *name);
static void setup_tests(void);
static void cleanup_tests(void);


int main(void)
{
  setup_tests();

  /*
    Whether or not to fork: forkin is wont kill the whole program if a test segfaults etc,
    but is difficult to debug a failing test.
   */
  enum fork_status fstat = CK_NOFORK;
  // enum fork_status fstat = CK_FORK;

  int number_failed;
  Suite *s;
  SRunner *sr;

  s = dct_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, fstat);

  // Only need to do this for ADDITIONAL suites. Otherwise it will run twice.


#if !defined(_USEMKL_)
  srunner_add_suite(sr, dct2_suite());
#endif
  srunner_add_suite(sr, cgsolve_suite());
  srunner_add_suite(sr, l1qc_newton_suite());
  srunner_add_suite(sr, vcl_math_suite());

  srunner_add_suite(sr, l1c_math_suite());
  srunner_add_suite(sr, l1c_common_suite());
  srunner_add_suite(sr, TV_suite());
  srunner_add_suite(sr, bregman_suite());

  /* Run the tests */
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  cleanup_tests();


  return(number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}



/*
  Setup path to data directory for all the tests.
 */
void setup_tests(void){

  /* We cannot free the result of getenv, but if the environmental variable did not exist,
   we should free the default we set.
  */
  char *srcdir;
  srcdir = getenv("srcdir");

  if (!srcdir){
    test_data_dir = fullfile(".", test_data_dir_name);
  }else{
    test_data_dir = fullfile(srcdir, test_data_dir_name);
  }
}

void cleanup_tests(void){
  free(test_data_dir);
}

/*
  Joins base_path with name as a file path and returns the result (malloced within).
  That is, it returns
  full_path = base_path/name.

  Memory for full_path will be allocated by this function.
  You are responsible for freeing full_path.
 */
char* fullfile(char *base_path, char *name){
  int len_base = strlen(base_path);
  int len_name = strlen(name);

  char *full_path = malloc( (len_base + len_name + 2)*sizeof(char));

  strncpy(full_path, base_path, len_base);
  full_path[len_base] = '/';
  strcpy(full_path + len_base + 1, name);

  return full_path;

}
