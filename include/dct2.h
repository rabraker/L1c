#ifndef _DCT2_
#define _DCT2_
#include "l1c_common.h"

extern int  dct2_setup(l1c_int N, l1c_int M, l1c_int Ny, l1c_int *pix_mask_idx);
extern void dct2_destroy();

extern void dct2_EMx(double *x_fftw, double *y);
extern void dct2_MtEty(double *y, double *x);
extern void dct2_idct(double *x_fftw);
extern void dct2_MtEt_EMx(double *x_fftw, double *z);
#endif
