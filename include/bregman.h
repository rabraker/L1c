#ifndef _BREGMAN_
#define _BREGMAN_
#include "l1c.h"

static inline double dsign(double a) { return a >= 0 ? 1.0 : -1.0; }

/** Pointer table to all the un-exposed functions in bregman.c so we can call
 from unit tests.
*/
typedef struct BregFuncs {
  void (*breg_shrink1)(l1c_int N, double* x, double* d, double gamma);

  void (*breg_mxpy_z)(l1c_int N, double* restrict x, double* restrict y, double* z);

  void (*breg_anis_jacobi)(
      int n, int m, double* uk_1, double* uk, double* rhs, double* D, double lambda);

  void (*breg_anis_guass_seidel)(int n, int m, double* u, double* rhs, double mu, double lambda);

  void (*breg_anis_rhs)(l1c_int n,
                        l1c_int m,
                        double* f,
                        double* dx,
                        double* bx,
                        double* dy,
                        double* by,
                        double* rhs,
                        double mu,
                        double lambda,
                        double* dwork1,
                        double* dwork2);

  void (*hess_inv_diag)(l1c_int n, l1c_int m, double mu, double lambda, double* D);
} BregFuncs;

BregFuncs breg_get_functions();

// static void breg_minx_y_z(l1c_int N, double * restrict x, double * restrict y, double *z);

int breg_anistropic_TV(l1c_int n,
                       l1c_int m,
                       double* uk,
                       double* f,
                       double mu,
                       double tol,
                       int max_iter,
                       int max_jac_iter);

// void breg_rhs(l1c_int n, l1c_int m, double *f, double *dx, double *bx, double *dy, double *by,
// double *rhs, double mu, double lambda, double *dwork1, double *dwork2);

// void breg_hess_eval(l1c_int N, double *x, double *y, void *hess_data);

// void hess_inv_diag(l1c_int n, l1c_int m, double mu, double lambda, double *D);
#endif
