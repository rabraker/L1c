#include "config.h"
#include <stdbool.h>
#include <stdint.h>

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


/**
 * Allocate a **zeroed** array of doubles aligned to a DALGIN byte boundary.
 * @param[in] N Number of requested elements. The returned pointer should be
 * freed with l1c_free_double().
 *
 * @return pointer to an array of doubles of length N.
 *
 * @see l1c_malloc_double()
 * @see l1c_malloc_double()
 */
double* l1c_calloc_double(int N){
  double *dptr;
  dptr = l1c_malloc_double(N);
  if (dptr){
    l1c_init_vec(N, dptr, 0);
  }
  return dptr;
}
/**
 * Allocate a 2D array **zeroed** array of doubles with nrow rows and ncol
 * columns. Each row is aligned to a DALGIN byte boundary.
 *
 * @param[in] nrow Number of requested rows.
 * @param[in] ncol Number of requested elements in each row..
 *
 * @return pointer to an array of doubles of length N.
 *
 * @see l1c_malloc_double_2D()
 * @see l1c_free_double_2D()
 */
double** l1c_calloc_double_2D(l1c_int nrow, l1c_int ncol){
  double **dptr = l1c_malloc_double_2D(nrow, ncol);
  if (dptr){
    for(int i=0; i<nrow; i++){
      l1c_init_vec(ncol, dptr[i], 0);
    }
  }
  return dptr;
}

/**
 * Allocate a 2D array of unitialized doubles with nrow rows and ncol
 * columns. Each row is aligned to a DALGIN byte boundary.
 *
 * @param[in] nrow Number of requested rows.
 * @param[in] ncol Number of requested elements in each row..
 *
 * @return pointer to an array of doubles of length N.
 *
 * @see l1c_calloc_double_2D()
 * @see l1c_free_double_2D()
 */
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

  return NULL;

}

/**
 * Free a 2D array allocated by l1c_malloc_double_2D() or l1c_calloc_double_2D().
 * The function will first free each row, then the pointer containing the rows.
 *
 * @param[in] nrow Number of rows in the 2D array to be freed.
 * @param[in] ddptr An array (of size nrow) of pointers.
 */
void l1c_free_double_2D(int nrow, double **ddptr){

  /* If ddptr is NULL, we cant free the inner arrays.*/
  if (ddptr){
    for (int k=0; k<nrow; k++){
      l1c_free_double(ddptr[k]);
    }
  }
  free(ddptr);
}

/*
  Check if the double pointer p is aligned to DALIGN boundary.
  If yes, return true, else return false.
 */
bool is_aligned(double *p){
  uintptr_t ptr_as_int = (uintptr_t)p;

  if (0 == (ptr_as_int & (DALIGN-1)) ){
    return true;
  }
  return false;
}


/*
  Will return an interger n such that
  p+n is aligned to a DALIGN byte boundary, ie,

  ((p+n)/DALIGN)*DALIGN == p+n

 */
int next_daligned_offset(void *p){
  if (is_aligned(p))
    return 0;

  uintptr_t ptr_as_int = (uintptr_t)p;
  uintptr_t offset = DALIGN - (ptr_as_int & (DALIGN - 1));

  return (int)(offset)/sizeof(double);
}
