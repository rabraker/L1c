#ifndef _L1C_TRANSFORMS_
#define _L1C_TRANSFORMS_
#include "config.h"
#include "l1c_common.h"


typedef struct L1cAxFuns {
  void(*Ax)(double *x, double *y);
  void(*Aty)(double *y, double *x);
  void(*AtAx)(double *x, double *z);

  /* Optional. Used by dct and dct2.*/
  void(*M)(double *x);
  void(*Mt)(double *y);
  /* Currently unused.*/

  void(*E)(double *x);
  void(*Et)(double *y);

  /* All implementations must define destroy*/
  void(*destroy)(void);

  /* Reserved for use.*/
  void *data;
}L1cAxFuns;


int dct2_setup(l1c_int Nx, l1c_int Mx, l1c_int Ny, l1c_int *pix_mask_idx,  L1cAxFuns *ax_funs);

int dct1_setup(l1c_int Nx, l1c_int Ny, l1c_int *pix_mask_idx, L1cAxFuns *ax_funs);

int setup_matrix_transforms(l1c_int n, l1c_int m, double *A, L1cAxFuns *ax_funs);

int l1c_setup_dct_transforms(l1c_int Nx, l1c_int Mx, l1c_int Ny,
                             l1c_int *pix_idx, L1cAxFuns *ax_funs);
#endif
