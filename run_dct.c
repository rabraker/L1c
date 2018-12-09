#include "dct.h"

void main(){
  int i;
  int N = 512*512;
  double To = 0.25;
  double Ts = 1/50.0;
  double omega = 2.0*PI/To;

  double *x, *y, *x_dct_x;
  x = fftw_alloc_real(N);
  // y = fftw_alloc_real(N);


  // int pix_mask_len = 10;
  // int pix_mask_idx[] = {1, 10, 15, 20, 25, 30, 35, 40, 45, 49};
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
    x[i] = sin(omega * (double)i * Ts);
    x_dct_x[i] = sin(omega * (double)i * Ts);
    // y[i] = sin(omega * (double)i * Ts);
    // printf("val = %f\n", val);
  }

  clock_t begin = clock();
  // Compute (M^T * E^T ) * E * M *x
  // for (i=1; i<2; i++){
    y = dct_EMx(x_dct_x);
    x_dct_x = dct_MtEty(y);
  // }
  clock_t end = clock();

  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("total time: %f\n", time_spent);


  for(i=0; i< min(50, pix_mask_len); i++){
    // printf("x[%d]=%f,    y[%d]=%f\n", i, x[i],i, y[i]);
    printf("x[%d]=%f,    y[%d]=%f\n", i, x_dct_x[i],i, x[i]);
  }

  dct_destroy();
  // free(pix_mask_idx);
  free(x);
}
