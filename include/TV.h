#ifndef _L1C_TV_
#define _L1C_TV_
#include "config.h"
#include "l1c.h"

void l1c_Dx(l1c_int n, l1c_int m, double alpha, double *X, double *Dx);

void l1c_DxT(l1c_int n, l1c_int m, double alpha, double *A, double *dxt);

void l1c_DxTDx(l1c_int n, l1c_int m, double alpha, double *A, double *dxtdx);

void l1c_Dy(l1c_int n, l1c_int m, double alpha, double *X, double *Dy);

void l1c_DyT(l1c_int n, l1c_int m, double alpha, double *A, double *dyt);

void l1c_DyTDy(l1c_int n, l1c_int m, double alpha, double *A, double *dytdy);


#endif
