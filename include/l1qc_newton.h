/* Function defs for l1qc_newton.c

 */

#ifndef _L1QC_NEWTON_
#define _L1QC_NEWTON_
#include "config.h"

#include "cgsolve.h"
#include "l1c_common.h"
#include "l1c_transforms.h"


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

typedef struct NewtParams{
  double epsilon;
  double tau;
  double mu;
  double newton_tol;
  l1c_int newton_max_iter;
  l1c_int lbiter;
  double lbtol;
  double l1_tol;
  l1c_int verbose;
  CgParams cg_params;
  int warm_start_cg; /** 0: no WS; 1 use dx; 2: use step*dx */

}NewtParams;

typedef struct LBResult{
  double l1;
  int    total_newton_iter;
  int    total_cg_iter;
  int    status;

}LBResult;


double find_max_step(l1c_int N, GradData gd, double *fu1,
                     double *fu2, int M, double *r, double *DWORK,
                     double epsilon, L1cAxFuns Ax_funs);

LSStat line_search(l1c_int N, l1c_int M, double *x, double *u, double *r, double *b,
                   double *fu1, double *fu2, GradData gd,
                   LSParams ls_params, double **DWORK5, double *fe, double *f, L1cAxFuns Ax_funs);

void l1qc_hess_grad(l1c_int N, double *fu1, double *fu2, double *sigx, double *atr,
                    double fe,  double tau, GradData gd);
int l1qc_descent_dir(l1c_int N, double *fu1, double *fu2, double *r, double fe, double tau,
                    GradData gd, double **Dwork7, CgParams cg_params, CgResults *cg_result,
                    L1cAxFuns Ax_funs);

void H11pfun(l1c_int N, double *z, double *y,  void *hess_data_in);

int newton_init(l1c_int N, double *x, double *u,  NewtParams *params);


/* Evalutes the value function */
void f_eval(l1c_int N, double *x, double *u, l1c_int M, double *r, double *b,
            double tau, double epsilon, double *fu1, double *fu2,
            double *fe, double *f, L1cAxFuns Ax_funs);
LBResult l1qc_newton(l1c_int N, double *x, l1c_int M, double *b,
                            NewtParams params, L1cAxFuns Ax_funs);
#endif
