#ifndef __L1C__
#define __L1C__
// #include "config.h"
typedef int l1c_int;

typedef struct LBResult{
  double l1;
  int    total_newton_iter;
  int    total_cg_iter;
  int    status;

}LBResult;

typedef struct L1qcDctOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  double lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cgtol;
  int cgmaxiter;
  int warm_start_cg;
}L1qcDctOpts;


int l1qc_dct(int Nrow, int Ncol, double *x_out, int M, double *b, l1c_int *pix_idx,
              L1qcDctOpts opts, LBResult *lb_res);

int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
                       double mu, double tol, int max_iter, int max_jac_iter);
#endif
