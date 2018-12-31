#ifndef _DCT_
#define _DCT_

extern void dct_setup(int Nx, int Ny, int *pix_mask_idx);
extern void dct_destroy();

extern void dct_EMx_new(double *x_fftw, double *y);
extern void dct_MtEty(double *y, double *x);
extern double* idct_plain(double *x_fftw);
extern void dct_MtEt_EMx_new(double *x_fftw, double *z);
#endif
