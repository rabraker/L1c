#ifndef _L1C_LOGGING_
#define _L1C_LOGGING_ 1
#include <stdio.h>

typedef int l1c_printf_t(const char* format, ...);

extern l1c_printf_t* l1c_printf;

l1c_printf_t* l1c_replace_printf(l1c_printf_t* new_printf);

#endif
