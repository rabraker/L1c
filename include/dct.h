#ifndef _DCT_
#define _DCT_

extern void dct_setup(int Nx, int Ny, int *pix_mask_idx);
extern void dct_destroy();

extern double* dct_EMx_new(double *x_fftw);
extern double* dct_EMx();
extern double* dct_MtEty(double *y);
extern void dct_load_x(double *x);
extern double* dct_x_ptr();

extern double* dct_MtEt_EMx_new(double *x);
#endif
