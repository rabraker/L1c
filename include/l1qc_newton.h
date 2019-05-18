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
  double *Dwork_1m;
  void(*AtAx)(double *x, double *z);
}Hess_data;

typedef struct LineSearchParams {
  double s;
  double alpha;
  double beta;
} LSParams;

typedef struct LSStat {
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

typedef struct l1c_l1qcProb {
  /*Static problem data.*/
  l1c_int m;
  l1c_int n;
  double *b;
  double tau;
  double epsilon;
  /*Internal, dynamic data*/
  double *r;
  double *fu1;
  double *fu2;
  double fe_val;
  double f_val;
  /* gradient and hessian data*/
  double *w1p;
  double *dx;
  double *du;
  double *sig11;
  double *sig12;
  double *ntgu;
  double *gradf;
  double **DWORK7;

  l1c_AxFuns ax_funs;
}l1c_l1qcProb;


/* These could be static. But that makes them hard to test. They are complicated enough I think
 they should be tested on their own. I.e., I dont want to *just* test l1qc_newton()
*/
double _l1c_l1qc_find_max_step(l1c_l1qcProb *Prb);


void _l1c_l1qc_hess_grad(l1c_l1qcProb *l1qc_prob, double *sigx, double *atr);

int _l1c_l1qc_descent_dir(l1c_l1qcProb *l1qc_prob, l1c_CgParams cg_params,
                          l1c_CgResults *cg_result);

void _l1c_l1qc_H11pfun(l1c_int m, double *z, double *y,  void *hess_data_in);


int _l1c_l1qc_check_feasible_start(l1c_l1qcProb *Prbm, double *x);

/* Evalutes the value function */
double _l1c_l1qc_f_eval(void *prob_data, double *x, double *u);

int _l1c_l1qc_newton_init(l1c_int m, double *x, double *u,  l1c_L1qcOpts *params);

#endif
