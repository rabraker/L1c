#ifndef _L1QC_COMMON_
#define _L1QC_COMMON_
#include "config.h"

#define DALIGN  64
#define ALIGNMENT_DOUBLE DALIGN

#include <stdio.h>


typedef int l1c_int;

#if defined(_USEMKL_)
   #include "mkl.h"
#else
   /* If we dont have MKL, use standard cblas.h. */
   #include "cblas.h"
#endif // _USE_MKL


static inline double min(double a,double b){
    return a < b ? a : b;
}

static inline double max(double a, double b){
    return a > b ? a : b;
}

/*
#define max(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a > _b ? _a : _b; })
*/

/*
  #define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
  __typeof__ (b) _b = (b);                    \
  _a < _b ? _a : _b; })
*/

#define L1C_INFEASIBLE_START  (1U << 1)
#define L1C_OUT_OF_MEMORY     (1U << 3)
#define L1C_DCT_INIT_FAILURE  (1U << 5)
#define L1C_FILE_READ_FAILURE (1U << 7)
#define L1C_CGSOLVE_FAILURE   (1U << 7)


#endif
