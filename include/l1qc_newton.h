/* Function defs for l1qc_newton.c

 */

#ifndef _L1QC_NEWTON_
#define _L1QC_NEWTON_


typedef struct Hess_data_ {
  double one_by_fe;
  double one_by_fe_sqrd;
  double *atr;
  double *sigx;

}Hess_data;

double sum_vec(int N, double *x);
void log_vec(int N, double alpha, double *x, double *logx);


void H11pfun(int N, double *z, double *y,  void *hess_data_in);
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
