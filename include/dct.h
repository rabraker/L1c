#ifndef _DCT_
#define _DCT_
#include "l1c_common.h"

extern int  dct_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx);
extern void dct_destroy();

extern void dct_EMx(double *x_fftw, double *y);
extern void dct_MtEty(double *y, double *x);
extern void dct_idct(double *x_fftw);
extern void dct_MtEt_EMx(double *x_fftw, double *z);
#endif
