#ifndef _L1C_MATH_
#define _L1C_MATH_
#include "config.h"
#include "l1c.h"

static inline double imin(l1c_int a, l1c_int b) { return a < b ? a : b; }

static inline double imax(l1c_int a, l1c_int b) { return a > b ? a : b; }

static inline double min(double a, double b) { return a < b ? a : b; }

static inline double max(double a, double b) { return a > b ? a : b; }

void __l1c_daxpby(l1c_int N, double alpha, double* restrict x, double beta, double* restrict y);

void l1c_init_vec(l1c_int N, double* x, double alpha);

void l1c_dxmuly_z(l1c_int N, double* x, double* y, double* z);

void l1c_daxpy_z(l1c_int N, double alpha, double* x, double* y, double* z);

void l1c_daxpby_z(l1c_int N, double alpha, double* x, double beta, double* y, double* z);

double l1c_dnorm1(l1c_int N, double* x);

double l1c_dsum(l1c_int N, double* x);

double l1c_dlogsum(l1c_int N, double alpha, double* x);

double l1c_dnrm2_rel_err(l1c_int N, double* restrict x, double* restrict y);

double l1c_dnrm2_err(l1c_int N, double* restrict x, double* restrict y);

double l1c_max_vec(l1c_int n, double* x);
void l1c_abs_vec(l1c_int n, double* x, double* y);

#if defined(HAVE_CBLAS_DAXPBY)
#define l1c_daxpby(N, alpha, x, beta, y) cblas_daxpby((N), (alpha), (x), 1, (beta), (y), 1)
#elif defined(HAVE_ATLAS_DAXPBY)
#define l1c_daxpby(N, alpha, x, beta, y) catlas_daxpby((N), (alpha), (x), 1, (beta), (y), 1)
#else
#define l1c_daxpby(N, alpha, x, beta, y) __l1c_daxpby((N), (alpha), (x), (beta), (y))
#endif

#endif
