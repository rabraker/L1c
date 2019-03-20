#ifndef __L1C_MATH__
#define __L1C_MATH__
#include "config.h"
#include "l1c_common.h"


void __l1c_daxpby(l1c_int N,
                  double alpha,
                  double * restrict x,
                  double beta,
                  double * restrict y);

void l1c_init_vec(l1c_int N, double *x, double alpha);

void l1c_dxmuly_z(l1c_int N,
                  double *x,
                  double *y,
                  double *z);

void l1c_daxpy_z(l1c_int N,
                 double alpha,
                 double *x,
                 double *y,
                 double *z);

void axpby_z(l1c_int N,
             double alpha,
             double *x,
             double beta,
             double *y,
             double *z);

double l1c_dnorm1(l1c_int N,
                  double *x);

double l1c_dsum(l1c_int N,
                double *x);

double l1c_dlogsum(l1c_int N,
                   double alpha,
                   double *x);


#if defined(HAVE_CBLAS_DAXPBY)
#define l1c_daxpby(N, alpha, x, beta, y) cblas_daxpby((N), (alpha), (x), 1, (beta), (y), 1)
#elif defined(HAVE_ATLAS_DAXPBY)
#define l1c_daxpby(N, alpha, x, beta, y) catlas_daxpby((N), (alpha), (x), 1, (beta), (y), 1)
#else
#define l1c_daxpby(N, alpha, x, beta, y) __l1c_daxpby((N), (alpha), (x), (beta), (y))
#endif

#endif
