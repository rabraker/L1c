#ifndef _DCTMKL_
#define _DCTMKL_
#include "l1c_common.h"

extern int dctmkl_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx);
extern void dctmkl_destroy();

extern void dctmkl_idct(double * restrict x, double * restrict y);

extern void dctmkl_EMx_new(double * restrict x_fftw, double * restrict y);
extern void dctmkl_MtEty(double *y, double *x);
extern void dctmkl_MtEt_EMx_new(double *x_fftw, double *z);
#endif
