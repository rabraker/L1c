#ifndef _L1C_LINESEARCH_
#define _L1C_LINESEARCH_

#include "l1c.h"
#define MAX_LINESEARCH_ITER 32

LSStat l1c_linesearch_xu(l1c_int m, double *x, double *u, double *dx, double *du,
                             double *fx0, double g_dot_dxu, void *prob_data,
                             double (*feval)(void *data, double *x, double *u),
                             LSParams ls_params, double **DWORK5);

LSStat l1c_linesearch(l1c_int n, double *x, double *dx,
                      double *fx0, double g_dot_dx, void *obj_data,
                      double (*objective_fun)(void *data, double *x),
                      LSParams ls_params, double *DWORK);
#endif
