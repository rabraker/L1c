#ifndef __VCL_MATH__
#define __VCL_MATH__
#include "l1c_common.h"

#ifdef __cplusplus
extern "C" {
#endif
  double vcl_logsum(const l1c_int N, const double alpha, const double *x);

  double vcl_sum(const l1c_int N, const double *x);
#ifdef __cplusplus
}
#endif

#endif
