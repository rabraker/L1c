#ifndef __L1C_NESTA__
#define __L1C_NESTA__

#define L1C_NESTA_NMEAN 10


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
  /** (final) Accuracy*/
  double mu;
  /** (final) Termination tolerance.*/
  double tol;
  /** Accuracy at continuation iteration j*/
  double mu_j;
  /** (final) Termination tolerance at continuation iteration j*/
  double tol_j;

  double L_by_mu;

  int n_continue;

  unsigned flags;

} l1c_NestaProb;


struct l1c_fmean_fifo {
  double *f_vals;
  int n_total;
  double *next; /*Pointer to the next location to write to.  */
};



void l1c_free_nesta_problem(l1c_NestaProb *NP);

l1c_NestaProb* _l1c_NestaProb_new(l1c_int n, l1c_int m);

void l1c_nesta_project(l1c_NestaProb *NP, double *xx, double *g, double *vk);

void l1c_nesta_feval(l1c_NestaProb *NP);

int l1c_nesta_setup(l1c_NestaProb *NP, double *beta_mu, double *beta_tol,
                    int n_continue, double *b, l1c_AxFuns ax_funs, double sigma,
                    double mu, double tol, double L, unsigned flags);

struct l1c_fmean_fifo _l1c_new_fmean_fifo(void);
void _l1c_push_fmeans_fifo(struct l1c_fmean_fifo *fifo, double fval);
double _l1c_mean_fmean_fifo(struct l1c_fmean_fifo *fifo);

int l1c_nesta(l1c_int m, double *xk, double mu, l1c_int n, double *b,
              l1c_AxFuns ax_funs, double sigma);


#define L1C_SYNTHESIS (1U << 0)
#define L1C_ANALYSIS (1U << 1)
#define L1C_WITH_DV (1U << 2)
#define L1C_WITH_DH (1U << 3)


#endif
