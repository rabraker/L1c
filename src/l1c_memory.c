#include "config.h"
#include "l1c_common.h"

#if defined(_HAVE_POSIX_MEMALIGN_)
// #define _POSIX_C_SOURCE  200112L
#include <stdlib.h>

double* malloc_double(int N){
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

void free_double(double *x){
  free(x);
}


#elif defined(_HAVE_MM_MALLOC_)
#include <xmmintrin.h>

double* malloc_double(int N){
  void *dptr;
  dptr = _mm_malloc((size_t)(N) * sizeof(double), DALIGN);
  return (double*) dptr;
}

void free_double(double *x){
  _mm_free(x);
}
#endif //_HAVE_POSIX_MEMALIGN_
