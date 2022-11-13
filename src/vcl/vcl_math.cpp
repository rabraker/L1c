#include <math.h>
#include <stdio.h>

#define VCL_NAMESPACE vcl
#include "vcl/vectorclass.h"
#include "vcl/vectormath_exp.h"

using namespace vcl;

#define VECTORSIZE 4

extern "C" double vcl_logsum(const int N, const double alpha, const double* x) {

  const int regularpart = N & (~(VECTORSIZE - 1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i = 0;
  double total = 0.0;
  Vec4d xvec, sum_vec(0);
  // if (N>4){
  for (i = 0; i < regularpart; i += VECTORSIZE) {
    xvec.load_a(x + i);
    sum_vec += vcl::log(xvec * alpha);
  }
  // }
  for (; i < N; i++) {
    total += log(x[i] * alpha);
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}

extern "C" double vcl_sum(const int N, const double* x) {
  const int regularpart = N & (~(VECTORSIZE - 1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i = 0;
  double total = 0.0;
  Vec4d xvec, sum_vec(0.0);
  for (i = 0; i < regularpart; i += VECTORSIZE) {
    xvec.load_a(x + i);
    sum_vec += xvec;
  }

  for (; i < N; i++) {
    total += x[i];
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}

/**
   Computes
   z = alpha*x
 */
extern "C" void vcl_dxMy_pz(const int N, const double* x, const double* y, double* z) {
  const int regularpart = N & (~(VECTORSIZE - 1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i = 0;
  Vec4d xvec, yvec, zvec;
  // #pragma omp parallel for private(xvec, yvec, zvec, i)
  for (i = 0; i < regularpart; i += VECTORSIZE) {
    xvec.load_a(x + i);
    yvec.load_a(y + i);
    zvec.load_a(z + i);

    zvec = xvec * yvec + zvec;
    zvec.store(z + i);
  }

  for (; i < N; i++) {
    z[i] = x[i] * y[i] + z[i];
  }
}
