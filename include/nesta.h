#ifndef _L1C_NESTA_
#define _L1C_NESTA_

#include "l1c.h"

#define L1C_NESTA_NMEAN 10
#define L1C_NESTA_MAX_INNER_ITER 10000

typedef struct _l1c_NestaProb l1c_NestaProb;

/**
 * Problem instance for Nesta.
 */
struct _l1c_NestaProb {
  /** Width of `R`, height of `W`. */
  l1c_int m;
  /** height of `R` */
  l1c_int n;
  /** Width of `W`*/
  l1c_int p;
  /** Value of (smoothed) functional at current iteration*/
  double fx;
  /** Prox center at current continuation iteration. */
  double* xo;
  /** Iterate.*/
  double* xk;
  /** Iterate.*/
  double* yk;
  /** Iterate.*/
  double* zk;
  /** We only want to compute A^T(b) once, so store it.*/
  double* Atb;
  /** Gradient of smoothed functional at current iteration..*/
  double* gradf;
  /** Weighted sum of gradient history
   * \f$\sum_i \alpha_i \nabla f_{\mu}(x_i) \f$
   * See, e.g., (3.8)
   */
  double* gradf_sum;
  /** Measurements.*/
  double* b;
  /** Work array of size(max(m, p)).*/
  double* dwork1;
  /** Work array of size(max(m, p)).*/
  double* dwork2;
  /** Transforms structure.*/
  l1c_AxFuns ax_funs;
  /** Noise level, \f$ ||Ax-b||\leq \sigma \f$ */
  double sigma;
  /** (final) Accuracy*/
  double mu;
  /** (final) Termination tolerance.*/
  double tol;
  /** Accuracy at continuation iteration j*/
  double mu_j;
  /** (final) Termination tolerance at continuation iteration j*/
  double tol_j;
  /** L, is the lipschitz constant of W (ie U), ie \f$||W||_2 \f$ */
  double L;
  /** Number of continuation iterations*/
  int n_continue;
  /** One of `L1C_ANALYSIS` or `L1C_SYNTHESIS`*/
  unsigned flags;
};

struct l1c_fmean_fifo {
  double* f_vals;
  int n_total;
  double* next; /*Pointer to the next location to write to.  */
};

void l1c_free_nesta_problem(l1c_NestaProb* NP);

l1c_NestaProb* _l1c_NestaProb_new(l1c_AxFuns ax_funs);

void l1c_nesta_project(l1c_NestaProb* NP, double* xx, double* g, double* vk);

void l1c_nesta_feval(l1c_NestaProb* NP);

int l1c_nesta_setup(l1c_NestaProb* NP,
                    double* beta_mu,
                    double* beta_tol,
                    double* b,
                    l1c_AxFuns ax_funs,
                    l1c_NestaOpts* opts);

struct l1c_fmean_fifo _l1c_new_fmean_fifo(void);
void _l1c_push_fmeans_fifo(struct l1c_fmean_fifo* fifo, double fval);
double _l1c_mean_fmean_fifo(struct l1c_fmean_fifo* fifo);

int l1c_nesta(l1c_int m,
              double* xk,
              l1c_int n,
              double* b,
              l1c_AxFuns ax_funs,
              l1c_NestaOpts opts);

#endif
