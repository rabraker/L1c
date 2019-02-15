/* Function defs for l1qc_newton.c

 */

#ifndef _L1QC_NEWTON_
#define _L1QC_NEWTON_
#include "config.h"

#include "cgsolve.h"
#include "l1c_common.h"

typedef struct AxFuns {
  void(*Ax)(double *x, double *y);
  void(*Aty)(double *y, double *x);
  void(*AtAx)(double *x, double *z);
}AxFuns;

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
  double *Adx;


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
                     double *fu2, int M, double *r, double epsilon);

LSStat line_search(l1c_int N, l1c_int M, double *x, double *u, double *r, double *b,
                   double *fu1, double *fu2, GradData gd,
                   LSParams ls_params, double *DWORK_5N, double *fe, double *f, AxFuns Ax_funs);

void get_gradient(l1c_int N, double *fu1, double *fu2, double *sigx, double *atr,
                  double fe,  double tau, GradData gd);
int compute_descent(l1c_int N, double *fu1, double *fu2, double *atr, double fe, double tau,
                    GradData gd, double *Dwork_6N, CgParams cg_params, CgResults *cg_result,
                    AxFuns AxFuns);

void H11pfun(l1c_int N, double *z, double *y,  void *hess_data_in);

int newton_init(l1c_int N, double *x, double *u,  NewtParams *params);


/* Evalutes the value function */
void f_eval(l1c_int N, double *x, double *u, l1c_int M, double *r, double tau, double epsilon,
                   double *fu1, double *fu2, double *fe, double *f);
LBResult l1qc_newton(l1c_int N, double *x, l1c_int M, double *b,
                            NewtParams params, AxFuns Ax_funs);
#endif
