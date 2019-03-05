#ifndef __BREGMAN___
#define __BREGMAN__
#include "l1c_common.h"

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
#endif
