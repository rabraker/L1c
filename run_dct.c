#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

#include "dct.h"
#include "dct_mkl.h"
#include "l1c_common.h"

#include "fftw/fftw3_mkl.h"

#define PI 3.141592653589793
#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })


int main(){
  int i;
  // int N = 512*512;
  int N = 16;
  // double To = 0.25;
  // double Ts = 1/50.0;
  // double omega = 2.0*PI/To;

  double *x, *y_fftw, *y_mkl;
  // double eta;
  x = malloc_double(N);
  y_fftw = malloc_double(N);
  y_mkl = malloc_double(N);

  int N_trials = 100;
  // int pix_mask_len = 23617;
  int pix_mask_len = N;
  int *pix_mask_idx;
  pix_mask_idx = calloc(pix_mask_len, sizeof(int));

  for(i=0; i<pix_mask_len; i++){
    pix_mask_idx[i] = i;
  }

  /* ----------------------------- */
  dct_setup(N, pix_mask_len, pix_mask_idx);
  dctmkl_setup(N, pix_mask_len, pix_mask_idx);


  for (i=0; i<N; i++){
    // eta = (double) (rand()/32767);
    // x[i] = sin(omega * (double)i * Ts) + eta;
    x[i] = (double)i;
  }

  dctmkl_EMx_new(x, y_mkl);
  dct_EMx_new(x, y_fftw);

  clock_t begin_mkl = clock();
  // Compute (M^T * E^T ) * E * M *x
  for (i=1; i<=N_trials; i++){
    dctmkl_EMx_new(x, y_mkl);
  }
  clock_t end_mkl = clock();

  clock_t begin_fftw = clock();
  // Compute (M^T * E^T ) * E * M *x
  for (i=1; i<=N_trials; i++){
    dct_EMx_new(x, y_fftw);
  }
  clock_t end_fftw = clock();


  double time_spent_fftw = (double)(end_fftw - begin_fftw) / CLOCKS_PER_SEC;
  double time_spent_mkl = (double)(end_mkl - begin_mkl) / CLOCKS_PER_SEC;
  printf("fftw time: %f   mkl-time: %f\n", time_spent_fftw/(double) N_trials, time_spent_mkl/(double) N_trials);


  for(i=0; i< min(50, pix_mask_len); i++){
    // printf("x[%d]=%f,    y[%d]=%f\n", i, x[i],i, y[i]);
    printf("y_fftw[%d]=%f,    y_mkl[%d]=%f\n", i, y_fftw[i],i, y_mkl[i]);
  }

  dct_destroy();
  dctmkl_destroy();
  // free(pix_mask_idx);
  free_double(x);
}
