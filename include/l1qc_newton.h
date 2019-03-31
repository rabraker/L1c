/* Internal function and struct defs for l1qc_newton.c

   This file really only needs to get included by the testing code.

 */

#ifndef _L1QC_NEWTON_
#define _L1QC_NEWTON_
#include "config.h"

#include "l1c.h"


typedef struct Hess_data_ {
  double one_by_fe;
  double one_by_fe_sqrd;
  double *atr;
  double *sigx;
  double *Dwork_1N;
  void(*AtAx)(double *x, double *z);
}Hess_data;

typedef struct LineSearchParams {
  double s;
  double alpha;
  double beta;
  double epsilon;
  double tau;

} LSParams;

typedef struct LSStat {
  double flx;
  double flu;
  double flin;
  double step;
  int iter;
  int status;
}LSStat;

typedef struct GradData{
  double *w1p;
  double *dx;
  double *du;
  double *sig11;
  double *sig12;
  double *ntgu;
  double *gradf;
}GradData;

/* These could be static. But that makes them hard to test. They are complicated enough I think
 they should be tested on their own. I.e., I dont want to *just* test l1qc_newton()
*/
double _l1c_l1qc_find_max_step(l1c_int N, GradData gd, double *fu1,
                               double *fu2, int M, double *r, double *DWORK,
                               double epsilon, l1c_AxFuns Ax_funs);

LSStat _l1c_l1qc_line_search(l1c_int N, l1c_int M, double *x, double *u, double *r, double *b,
                             double *fu1, double *fu2, GradData gd, LSParams ls_params,
                             double **DWORK5, double *fe, double *f, l1c_AxFuns Ax_funs);

void _l1c_l1qc_hess_grad(l1c_int N, double *fu1, double *fu2, double *sigx, double *atr,
                         double fe,  double tau, GradData gd);

int _l1c_l1qc_descent_dir(l1c_int N, double *fu1, double *fu2, double *r, double fe, double tau,
                          GradData gd, double **Dwork7, l1c_CgParams cg_params, l1c_CgResults *cg_result,
                          l1c_AxFuns Ax_funs);

void _l1c_l1qc_H11pfun(l1c_int N, double *z, double *y,  void *hess_data_in);

/* Evalutes the value function */
void _l1c_l1qc_f_eval(l1c_int N, double *x, double *u, l1c_int M, double *r, double *b,
                      double tau, double epsilon, double *fu1, double *fu2, double *fe,
                      double *f, l1c_AxFuns Ax_funs);

int _l1c_l1qc_newton_init(l1c_int N, double *x, double *u,  l1c_L1qcOpts *params);

#endif
