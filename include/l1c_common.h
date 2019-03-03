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

#if defined(_HAVE_POSIX_MEMALIGN_)
   // #define _POSIX_C_SOURCE  200112L
   #include <stdlib.h>

   static inline double* malloc_double(int N){
     void *dptr;
     /*DALIGN must be multiple of sizeof(double) and power of two.
       This is satisfired for DALIGN=64 and sizeof(double)=8.
     */
     if(posix_memalign(&dptr, DALIGN, (size_t)(N) * sizeof(double))){
       return NULL; // We could check the value,
                 // it's different for out of Mem, vs unacceptable DALIGN.
     }else{
       return (double*) dptr;
     }
   }

   #define free_double(x) free(x)
#elif defined(_HAVE_MM_MALLOC_)
   #include <xmmintrin.h>
   #define malloc_double(N) (double*) _mm_malloc((size_t)(N) * sizeof(double), DALIGN)
   #define free_double(x) _mm_free(x)
#endif //_HAVE_POSIX_MEMALIGN_




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

#define L1C_INFEASIBLE_START  (1U << 1)
#define L1C_OUT_OF_MEMORY     (1U << 3)
#define L1C_DCT_INIT_FAILURE  (1U << 5)
#define L1C_FILE_READ_FAILURE (1U << 7)
#define L1C_CGSOLVE_FAILURE   (1U << 7)


#endif
