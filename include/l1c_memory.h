#ifndef _L1C_MEMORY_
#define _L1C_MEMORY_
#include "config.h"
#include "stdbool.h"


bool is_aligned(double *p);

int next_daligned_offset(void *p);



#endif
