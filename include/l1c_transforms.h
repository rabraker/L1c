#ifndef _L1C_TRANSFORMS_
#define _L1C_TRANSFORMS_
#include "config.h"
#include "l1c_common.h"

typedef struct L1cAxFuns {
  void(*Ax)(double *x, double *y);
  void(*Aty)(double *y, double *x);
  void(*AtAx)(double *x, double *z);

  /* Currently unused.*/
  void(*M)(double *x);
  void(*Mt)(double *y);
  void(*E)(double *x);
  void(*Et)(double *y);
  void(*destroy)(void);

  void *data;
}L1cAxFuns;



int setup_matrix_transforms(l1c_int n, l1c_int m, double *A, L1cAxFuns *ax_funs);
void destroy_matrix_transforms(void);
void Ax(double *x, double *y);
void Aty(double *y, double *x);
void AtAx(double *x_in, double *x_out);


#endif
