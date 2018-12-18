/* Function defs for l1qc_newton.c

 */

#ifndef _L1QC_NEWTON_
#define _L1QC_NEWTON_

#include "cgsolve.h"

typedef struct Hess_data_ {
  double one_by_fe;
  double one_by_fe_sqrd;
  double *atr;
  double *sigx;

}Hess_data;

typedef struct LineSearchParams {
  double s;
  double alpha;
  double beta;
  double epsilon;
  double tau;

} LSParams;

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
  double newton_max_iter;
  double lbiter;
  double lbtol;
  int verbose;
  CgParams cg_params;

}NewtParams;

double sum_abs_vec(int N, double *x);

double sum_vec(int N, double *x);

void log_vec(int N, double alpha, double *x, double *logx);

double logsum(int N, double *x, double alpha);

double find_max_step(int N, GradData gd, double *fu1,
                     double *fu2, int M, double *r, double epsilon, int *Iwork_2N);

int line_search(int N, int M, double *x, double *u, double *r, double *fu1, double *fu2, GradData gd,
                LSParams ls_params, double *DWORK_5N, double *fe, double *f);

void get_gradient(int N, double *fu1, double *fu2, double *sigx, double *atr,
                  double fe,  double tau, GradData gd);
int compute_descent(int N, double *fu1, double *fu2, double *r, double fe, double tau,
                    GradData gd, double *Dwork_5N, CgParams cg_params, CgResults *cg_result);

void H11pfun(int N, double *z, double *y,  void *hess_data_in);

int newton_init(int N, double *x, double *u,  NewtParams *params,
                double *Dwork,  int M, double *b, int *pix_idx);
/*
int get_gradient(int N, double *fu1, double *fu2, double fe,  double tau, double *gradf);
*/

/* Evalutes the value function */
extern void f_eval(int N, double *r, double *x, double *u, double tau, double epsilon,
                   double *fu1, double *fu2, double *fe, double *f, double *Dwork);
// void l1qc_newton(int N, double *x, double u,
//                  void A(int N, int m, double *x, double * y),
//                  void At(int N, int m, double *y, double * x),
//                  double *b, double epsilon, double tau, double newtontol,
//                  int newtonmaxiter, double cgtol, int cgmaxiter, inter verbose)

#endif
