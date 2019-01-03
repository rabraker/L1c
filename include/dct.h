#ifndef _DCT_
#define _DCT_
#include "l1qc_common.h"

extern int  dct_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx);
extern void dct_destroy();

extern void dct_EMx_new(double *x_fftw, double *y);
extern void dct_MtEty(double *y, double *x);
extern double* idct_plain(double *x_fftw);
extern void dct_MtEt_EMx_new(double *x_fftw, double *z);
#endif
