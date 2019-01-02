#ifndef _DCTMKL_
#define _DCTMKL_

extern void dctmkl_setup(int Nx, int Ny, int *pix_mask_idx);
extern void dctmkl_destroy();

extern void dctmkl_idct(double * restrict x, double * restrict y);

extern void dctmkl_EMx_new(double * restrict x_fftw, double * restrict y);
extern void dctmkl_MtEty(double *y, double *x);
extern void dctmkl_MtEt_EMx_new(double *x_fftw, double *z);
#endif
