#include "config.h"

#include <math.h>

#include "cblas.h"
#include "l1c.h"

/**
   Compute the vector sum
   y = alpha*x + beta*y.

   This is a replacement for cblas_daxpby. If cblas_daxpby is availible, use it.
*/
void __l1c_daxpby(
    l1c_int N, double alpha, double* restrict x, double beta, double* restrict y) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
#pragma omp parallel for
  for (l1c_int i = 0; i < N; i++) {
    y_[i] = alpha * x_[i] + beta * y_[i];
  }
}

/**
   Initialize a vector x of length N to alpha in all entries.
*/
void l1c_init_vec(l1c_int N, double* restrict x, double alpha) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    x_[i] = alpha;
  }
}

/**
   Compute the vector elementwise product
   z = x.*y
*/
void l1c_dxmuly_z(l1c_int N,
                  double* restrict x,
                  double* restrict y,
                  double* restrict z) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  double* z_ = __builtin_assume_aligned(z, DALIGN);

#pragma omp parallel for
  for (l1c_int i = 0; i < N; i++) {
    z_[i] = x_[i] * y_[i];
  }
}

/**
   Computes vector addition

   z = alpha*x + beta*y. Similar to cblas_daxpy, but for the situation when the result
   should not overwrite y.

*/
void l1c_daxpby_z(l1c_int N,
                  double alpha,
                  double* restrict x,
                  double beta,
                  double* restrict y,
                  double* restrict z) {
  /* Computes z = a * x + y. Similary to cblas_axpy, but for when you don't want to
     overwrite y. This way, we avoid a call to cblas_dcopy().

     May be worth explicitly vectorizing (e.g., ispc??) this function.
  */
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  double* z_ = __builtin_assume_aligned(z, DALIGN);

  l1c_int i;
#pragma omp parallel for
  for (i = 0; i < N; i++) {
    z_[i] = alpha * x_[i] + beta * y_[i];
  }
}

/**
Computes vector addition

z = a*x + y. Similar to cblas_daxpy, but for the situation when the result should not
overwrite y.

*/
void l1c_daxpy_z(l1c_int N,
                 double alpha,
                 double* restrict x,
                 double* restrict y,
                 double* restrict z) {
  /* Computes z = a * x + y. Similary to cblas_axpy, but for when you don't want to
   overwrite y. This way, we avoid a call to cblas_dcopy().

   May be worth explicitly vectorizing (e.g., ispc??) this function.
  */
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  double* z_ = __builtin_assume_aligned(z, DALIGN);

  l1c_int i;
#pragma omp parallel for
  for (i = 0; i < N; i++) {
    z_[i] = alpha * x_[i] + y_[i];
  }
}

/**
   Computes the 1-norm of the vector x, i.e.,
   dnrm1 = sum_0^{n-1} |x[i]|
 */
double l1c_dnorm1(l1c_int N, double* x) { return cblas_dasum(N, x, 1); }

double l1c_dsum(l1c_int N, double* x) {
  l1c_int i = 0;
  double sum = 0.0;
  for (i = 0; i < N; i++) {
    sum = sum + x[i];
  }
  return sum;
}

double l1c_dlogsum(l1c_int N, double alpha, double* x) {
  /* Computes sum(log( alpha *x)) */
  l1c_int i = 0;
  double total = 0.0;
#pragma omp parallel for reduction(+ : total)
  for (i = 0; i < N; i++) {
    total += log(alpha * x[i]);
  }
  return total;
}

/**
   Computes the relative 2-norm of error the error between x and y, i.e.,
   nrm = ||x - y||_2/||y||_2
   Assumes x and y are aligned to a DALIGN (default: 64) byte boundary.
*/
double l1c_dnrm2_rel_err(l1c_int N, double* restrict x, double* restrict y) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  int i = 0;
  double ynrm = cblas_ddot(N, y, 1, y, 1);
  double nrm = 0.0, diff = 0.0;

  for (i = 0; i < N; i++) {
    diff = x_[i] - y_[i];
    nrm += diff * diff;
  }
  return sqrt(nrm / ynrm);
}

/**
 * Computes the maximum element of the array x.
 */
double l1c_max_vec(l1c_int n, double* x) {

  if (n <= 0) {
    return NAN;
  }

  double xmax = x[0];
  if (n == 1) {
    return xmax;
  }

  for (int k = 1; k < n; k++) {
    xmax = x[k] > xmax ? x[k] : xmax;
  }

  return xmax;
}

/**
 * Computes the absoulte value of x and stores it in y.
 */
void l1c_abs_vec(l1c_int n, double* x, double* y) {

  for (int k = 0; k < n; k++) {
    y[k] = fabs(x[k]);
  }
}

/**
   Computes the 2-norm of error the error between x and y, i.e.,
   nrm = ||x - y||_2
   Assumes x and y are aligned to a DALIGN (default: 64) byte boundary.
*/
double l1c_dnrm2_err(l1c_int N, double* restrict x, double* restrict y) {
  double* x_ = __builtin_assume_aligned(x, DALIGN);
  double* y_ = __builtin_assume_aligned(y, DALIGN);
  int i = 0;
  double nrm = 0.0, diff = 0.0;

  for (i = 0; i < N; i++) {
    diff = x_[i] - y_[i];
    nrm += diff * diff;
  }
  return sqrt(nrm);
}
