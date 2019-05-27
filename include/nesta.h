#ifndef __L1C_NESTA__
#define __L1C_NESTA__



typedef struct l1c_NestaProb {
  l1c_int m;
  l1c_int n;
  /** */
  double fx;
  double *xo;
  double *xk;
  double *yk;
  double *zk;
  double *Atb;
  double *gradf;
  double *gradf_sum;
  double *b;
  double *dwork1;
  double *dwork2;
  l1c_AxFuns ax_funs;

  /** Noise level, \f||Ax-b||\leq sigma \f*/
  double sigma;
  /** Accuracy*/
  double mu;
  /** Termination tolerance.*/
  double tol;
  double L_by_mu;
  unsigned flags;
} l1c_NestaProb;


void l1c_free_nesta_problem(l1c_NestaProb *NP);

l1c_NestaProb* l1c_init_nesta_problem(l1c_int n, l1c_int m, double *b, l1c_AxFuns ax_funs,
                                      double sigma, double mu, double tol, double L, unsigned flags);

void l1c_nesta_project(l1c_NestaProb *NP, double *xx, double *g, double *vk);

void l1c_nesta_feval(l1c_NestaProb *NP);

#define L1C_SYNTHESIS (1U << 0)
#define L1C_ANALYSIS (1U << 1)
#define L1C_WITH_DV (1U << 2)
#define L1C_WITH_DH (1U << 3)


#endif
