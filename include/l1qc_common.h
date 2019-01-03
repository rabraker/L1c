#ifndef _L1QC_COMMON_
#define _L1QC_COMMON_

#define ALIGNMENT_DOUBLE 64

#ifdef __MATLAB__
#include "mex.h"
#else
#include <stdio.h>
#endif


#ifdef _USEMKL_
#include "mkl.h"
#define free_double(x) mkl_free(x)
// parenthesis around N are crucial here.
#define malloc_double(N) (double*)mkl_malloc((size_t)(N) *sizeof(double), ALIGNMENT_DOUBLE)
#else
#include "cblas.h"

#define free_double(x) free(x)
#define malloc_double(N) calloc(N, sizeof(double))
#endif

#if defined(MKL_INT)
typedef MKL_INT l1c_int;
#elif defined(blasint)
typedef blasint l1c_int;
#else
typedef int l1c_int;
#endif

/*
#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })
*/

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

#endif
