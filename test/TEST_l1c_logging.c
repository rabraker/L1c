#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <check.h>
#include "l1c_logging.h"
#include "json_utils.h"

static char fname_tmp[] = "test_scratch.txt";

int l1c_printf_new(const char *format, ...) {
  FILE *fid = fopen(fname_tmp, "w+");
  if (!fid){
    return -1;
  }
  fprintf(fid, "NEW:");

  va_list args;
  va_start(args, format);
  int status = vfprintf(fid, format, args);
  va_end(args);

  fclose(fid);
  return status;
}



 /* Check that we can reset the old printf */
START_TEST(test_l1c_replace_printf) {
  char str_expected[] = "NEW:test format: 5, 3.14";
  char *str_actual;
  l1c_printf_t *old_printf = l1c_replace_printf(l1c_printf_new);

  l1c_printf("test format: %d, %.2f", 5, 3.14);
  load_file_as_text(fname_tmp, &str_actual);
  ck_assert_str_eq(str_actual, str_expected);
  free(str_actual);

  /*Now, we re-open the file, write a single character,
    (which will erase the old contents) */
  FILE *fid = fopen(fname_tmp, "w+");
  if (!fid) {
    ck_abort_msg("Failed to open file");
  }
  fprintf(fid, "1");
  fclose(fid);

  /*Now, reset l1c_printf to the old_printf (which should be the system printf)*/
  l1c_replace_printf(old_printf);
  /* Now, when we call l1c_printf, it should not go to the file.
   I dont know how to check that it goes to stdout. */
  l1c_printf("test format: %d, %.2f", 5, 3.14);
  load_file_as_text(fname_tmp, &str_actual);

  ck_assert_str_eq(str_actual, "1");

  free(str_actual);
}
END_TEST

Suite *l1c_logging_suite(void) {
  Suite *s;

  TCase *tc_logging;
  s = suite_create("logging_suite");

  tc_logging = tcase_create("logging");
  tcase_add_test(tc_logging, test_l1c_replace_printf);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_logging);

  return s;
}
