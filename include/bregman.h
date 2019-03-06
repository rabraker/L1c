#ifndef __BREGMAN___
#define __BREGMAN__
#include "l1c_common.h"


typedef struct BregTVHessData {
  l1c_int n;
  l1c_int m;
  double *dwork1;
  double *dwork2;
  double mu;
  double lambda;
}BregTVHessData;


static inline double dsign(double a){
  return a >= 0 ? 1.0 : -1.0;
}

/** Pointer table to all the un-exposed functions in bregman.c so we can call
 from unit tests.
*/
typedef struct BregFuncs{
  void(*breg_shrink1)(l1c_int N, double *x, double *d, double gamma);
  void(* breg_mxpy_z)(l1c_int N, double * restrict x, double * restrict y, double *z);
}BregFuncs;

BregFuncs breg_get_functions();



void breg_minx_y_z(l1c_int N, double * restrict x, double * restrict y, double *z);

int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
                       double lambda, double tol, l1c_int max_iter);

void breg_rhs(l1c_int n, l1c_int m, double *f, double *dx, double *bx, double *dy, double *by,
              double *rhs, double mu, double lambda, double *dwork1, double *dwork2);

void breg_hess_eval(l1c_int N, double *x, double *y, void *hess_data);
#endif
