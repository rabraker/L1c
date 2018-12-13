#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <math.h>

#include "dct.h"

#ifdef _USEMKL_
#include "fftw/fftw3_mkl.h"
#endif

#define PI 3.141592653589793
#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })

void main(){
  int i;
  int N = 512*512*20;
  double To = 0.25;
  double Ts = 1/50.0;
  double omega = 2.0*PI/To;

  double *x, *y, *x_dct_x;
  double eta;
  x = fftw_alloc_real(N);

  int N_trials = 10;
  int pix_mask_len = 23617;
  int *pix_mask_idx;
  pix_mask_idx = calloc(pix_mask_len, sizeof(int));
  for(i=0; i<pix_mask_len; i++){
    pix_mask_idx[i] = i;
  }

  /* ----------------------------- */
  dct_setup(N, pix_mask_len, pix_mask_idx);

  x_dct_x = dct_x_ptr();

  for (i=0; i<N; i++){
    eta = (double) (rand()/32767);
    x[i] = sin(omega * (double)i * Ts) + eta;
    x_dct_x[i] = sin(omega * (double)i * Ts) + eta;
    // y[i] = sin(omega * (double)i * Ts);
    // printf("val = %f\n", val);
  }

  clock_t begin = clock();
  // Compute (M^T * E^T ) * E * M *x
  for (i=1; i<N_trials; i++){
    y = idct_plain(x_dct_x);
    // x_dct_x = dct_MtEty(y);
  }
  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("total time: %f\n", time_spent);
  printf("time per DCT: %f\n", time_spent/(double) N_trials);

  for(i=0; i< min(50, pix_mask_len); i++){
    // printf("x[%d]=%f,    y[%d]=%f\n", i, x[i],i, y[i]);
    printf("x[%d]=%f,    y[%d]=%f\n", i, x_dct_x[i],i, y[i]);
  }

  dct_destroy();
  // free(pix_mask_idx);
  fftw_free(x);
}
