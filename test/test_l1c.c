
/*

  Entry point for the L1c test suite.

 */

#include "config.h"
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "TEST.h"

char* test_data_dir;
static char test_data_dir_name[] = "test_data";
static char _TEST_DATA_DIR[] = "TEST_DATA_DIR";

char* fullfile(char* base_path, char* name);
static void setup_tests(void);
static void cleanup_tests(void);

int main(void) {
  setup_tests();

  /*
    Whether or not to fork: forkin is wont kill the whole program if a test segfaults etc,
    but is difficult to debug a failing test.
   */
  enum fork_status fstat = CK_NOFORK;
  // enum fork_status fstat = CK_FORK;

  int number_failed;
  Suite* s;
  SRunner* sr;

  s = dct_suite();
  sr = srunner_create(s);
  srunner_set_fork_status(sr, fstat);

  // Only need to do this for ADDITIONAL suites. Otherwise it will run twice.

  srunner_add_suite(sr, dct2_suite());
  srunner_add_suite(sr, dctTV_suite());
  srunner_add_suite(sr, cgsolve_suite());
  srunner_add_suite(sr, l1qc_newton_suite());
  srunner_add_suite(sr, vcl_math_suite());
  srunner_add_suite(sr, l1c_logging_suite());
  srunner_add_suite(sr, l1c_math_suite());
  srunner_add_suite(sr, l1c_nesta_suite());
  srunner_add_suite(sr, TV_suite());
  srunner_add_suite(sr, bregman_suite());
  srunner_add_suite(sr, l1c_linesearch_suite());
  srunner_add_suite(sr, l1c_memory_suite());
  srunner_add_suite(sr, matrix_transform_suite());

  /* Run the tests */
  srunner_run_all(sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  cleanup_tests();

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

/*
  Setup path to data directory for all the tests.
 */
void setup_tests(void) {

  /* We cannot free the result of getenv, but if the environmental variable did not exist,
   we should free the default we set.
  */
  char* data_dir;
  data_dir = getenv(_TEST_DATA_DIR);

  if (!data_dir) {
    test_data_dir = fullfile(".", test_data_dir_name);
  } else {
    test_data_dir = fullfile(data_dir, "");
  }
}

void cleanup_tests(void) { free(test_data_dir); }

/*
  Joins base_path with name as a file path and returns the result (malloced within).
  That is, it returns
  full_path = base_path/name.

  Memory for full_path will be allocated by this function.
  You are responsible for freeing full_path.
 */
char* fullfile(char* base_path, char* name) {
  int len_base = strlen(base_path);
  int len_name = strlen(name);

  /*Add 2: 1 for \0 and one for '/', path separator. */
  char* full_path = malloc((len_base + len_name + 2) * sizeof(char));
  if (!full_path) {
    return NULL;
  }
  /*Includes \0*/
  strcpy(full_path, base_path);
  full_path[len_base] = '/';
  /*len_base+1 will overwrite \0 from base_path. */
  strcpy(full_path + len_base + 1, name);

  return full_path;
}
