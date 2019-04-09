#include "config.h"
#include "l1c.h"
#include "l1c_math.h"

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
  double *dptr;
  dptr = _mm_malloc((size_t)(N) * sizeof(double), DALIGN);
  return dptr;
}

void l1c_free_double(double *x){
  _mm_free(x);
}
#endif //_HAVE_POSIX_MEMALIGN_


double* l1c_calloc_double(int N){
  double *dptr;
  dptr = l1c_malloc_double(N);
  if (dptr){
    l1c_init_vec(N, dptr, 0);
  }
  return dptr;
}

double** l1c_calloc_double_2D(l1c_int nrow, l1c_int ncol){
  double **dptr = l1c_malloc_double_2D(nrow, ncol);
  if (dptr){
    for(int i=0; i<nrow; i++){
      l1c_init_vec(ncol, dptr[i], 0);
    }
  }
  return dptr;
}

/** This is primarily for the DWORK arrays. */
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


/**
   Free an 2D array allocated by l1c_malloc_double_2D
*/
void l1c_free_double_2D(int nrow, double **dwork){

  for (int k=0; k<nrow; k++){
    l1c_free_double(dwork[k]);
  }
  free(dwork);
}
