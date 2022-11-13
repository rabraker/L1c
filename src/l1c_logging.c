#include "l1c_logging.h"
#include <stdio.h>

l1c_printf_t* l1c_printf = printf;

l1c_printf_t* l1c_replace_printf(l1c_printf_t* new_printf) {
  l1c_printf_t* old_printf = l1c_printf;

  l1c_printf = new_printf;

  return old_printf;
}
