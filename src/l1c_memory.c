#include "config.h"
#include "l1c.h"

#if defined(_HAVE_POSIX_MEMALIGN_)
// #define _POSIX_C_SOURCE  200112L
#include <stdlib.h>

double* l1c_malloc_double(int N){
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

void l1c_free_double(double *x){
  free(x);
}


#elif defined(_HAVE_MM_MALLOC_)
#include <xmmintrin.h>

double* l1c_malloc_double(int N){
  void *dptr;
  dptr = _mm_malloc((size_t)(N) * sizeof(double), DALIGN);
  return (double*) dptr;
}

void l1c_free_double(double *x){
  _mm_free(x);
}
#endif //_HAVE_POSIX_MEMALIGN_


/* This is primarily for the DWORK arrays. */
double** l1c_malloc_double_2D(l1c_int nrow, l1c_int ncol){
  int k;
  double **dwork = malloc(nrow*sizeof(double*));
  if (!dwork){
    goto fail1;
  }

  for (k=0; k<nrow; k++){
    dwork[k] = NULL;
  }

  for (k=0; k<nrow; k++){
    dwork[k] = l1c_malloc_double(ncol);
    if (!dwork[k]){
      goto fail2;
    }
  }

  return dwork;

 fail2:
  for (k=0; k<nrow; k++){
    l1c_free_double(dwork[k]);
  }
 fail1:
  free(dwork);

  return dwork;

}

void l1c_free_double_2D(int nrow, double **dwork){

  for (int k=0; k<nrow; k++){
    l1c_free_double(dwork[k]);
  }
  free(dwork);
}
