#include <math.h>

#include "l1c_common.h"

/**
   Compute the vector elementwise product
   z = x.*y
*/
void l1c_dxmuly_z(l1c_int N, double * restrict x, double * restrict y, double * restrict z){
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);
  double *z_ = __builtin_assume_aligned(z, 64);

#pragma omp parallel for
  for (l1c_int i=0; i<N; i++){
    z_[i] = x_[i] * y_[i];
  }
}

/**
Computes vector addition

z = a*x + y. Similar to cblas_daxpy, but for the situation when the result should not overwrite y.

*/
void l1c_daxpy_z(l1c_int N, double alpha, double * restrict x, double * restrict y, double * restrict z){
  /* Computes z = a * x + y. Similary to cblas_axpy, but for when you don't want to overwrite y.
   This way, we avoid a call to cblas_dcopy().

   May be worth explicitly vectorizing (e.g., ispc??) this function.
  */
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);
  double *z_ = __builtin_assume_aligned(z, 64);

  l1c_int i;
#pragma omp parallel for
  for (i = 0; i<N; i++){
    z_[i] = alpha * x_[i] + y_[i];
  }
}

/**
   Computes the 1-norm of the vector x, i.e.,
   dnrm1 = sum_0^{n-1} |x[i]|
 */
double l1c_dnorm1(l1c_int N, double *x){
  return cblas_dasum(N, x, 1);
}


double l1c_dsum(l1c_int N, double *x){
  l1c_int i = 0;
  double sum = 0.0;
  for (i=0; i<N; i++){
    sum = sum + x[i];
  }
  return sum;
}


double l1c_dlogsum(l1c_int N,  double alpha, double *x) {
  /* Computes sum(log( alpha *x)) */
  l1c_int i = 0;
  double total = 0.0;
#pragma omp parallel for reduction(+:total)
  for (i=0; i<N; i++){
    total += log(alpha * x[i]);
  }
  return total;
}
