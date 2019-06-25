/** @file nesta.c
 *
 * This is an implementation of the l1-specialized version of Nesterovs algorithm
 * described in @cite becker_nesta_2011.
 *
 */
#include <stdlib.h>
#include <cblas.h>
#include <math.h>
#include "l1c.h"
#include "l1c_math.h"
#include "nesta.h"
#include "l1c_logging.h"

/*
  W is a possibly overcomplete dictionary, ie, W is m by p. Signal v is assumed sparse in
  domain of W. In general, x = W*v \in \mathbb{R}^m, and v\in\mathbb{R}^p. This gives rise
  to the synthesis and analysis prblems:

  1) \f$ min_{v}  ||v||_1 s.t. || b - R * W * v||_2 \leq \sigma \f$

  2) \f$ min_{x}  ||W^T x||_1 s.t. || b - R * x||_2 \leq \sigma \f$


  In this code, when flags|synthesis, we set \f$W=I_{m,m}\f$, and use ax_funs.A.

  In the case of analysys, we assume that ax_funs.E = R and ax_funs.M= W.
  The following four functions provide wrappers to facilite this logic.

  The goal is that, e.g., dct1 is fully compaitible with this scheme for both
  synthesis and analysis. It might be better to rename ax_funs.E to R, since the
  connotation is more general.
 */





/* The next three functions basically implement a circular buffer for storing the vector
 containing the last L1C_NESTA_NMEAN values of fx. This is probably over-engineered...
*/

/*
   Returns a new instance of l1c_fmean_fifo, with the vector of
   fvals initialized to zero.
 */
struct l1c_fmean_fifo _l1c_new_fmean_fifo(void){
  double *fvals = malloc(sizeof(double)*L1C_NESTA_NMEAN);

  for (int i=0; i<L1C_NESTA_NMEAN; i++){
    fvals[i] = 0.0;
  }

  struct l1c_fmean_fifo fifo = {.f_vals=fvals,
                     .n_total=0,
                     .next=fvals};
  return fifo;
}

/** @private Push a new value of fval into the fifo and increment fifo.n_total.
 *
 * @param [in,out] fifo An instance of l1c_fmean_fifo.
 * @param [in] fval the new value of the functional, to be pushed into the fifo.
 */
void _l1c_push_fmeans_fifo(struct l1c_fmean_fifo *fifo, double fval) {

  *fifo->next = fval;

  if (fifo->next + 1 > fifo->f_vals + L1C_NESTA_NMEAN-1){
    fifo->next = fifo->f_vals;
  }else{
    fifo->next++;
  }

  if(fifo->n_total < L1C_NESTA_NMEAN){
    fifo->n_total++;
  }
}

/** @private
 * Compute the mean of fifo.fval. This is a standard mean,
 * except it is normalized by fifo.n_total, not the length of the buffer.
 * This is eq (3.13) in the paper.
 *
 * @param [in] fifo
 */
double _l1c_mean_fmean_fifo(struct l1c_fmean_fifo *fifo) {

  double mean = 0;
  /* Always run the sum over the whole thing. Unused elements should be zero.*/
  for (int i=0; i < L1C_NESTA_NMEAN; i++){
    mean += fifo->f_vals[i];
  }

  return (mean / (double)fifo->n_total);
}


/** @private
   Initialize a l1c_NestaProblem struct. Memory will be allocated for all
   arrays.
 */
l1c_NestaProb* _l1c_NestaProb_new(l1c_AxFuns ax_funs){

  l1c_NestaProb *NP = malloc(sizeof(l1c_NestaProb));
  if (!NP){
    return NULL;
  }

  l1c_int n = ax_funs.n;
  l1c_int m = ax_funs.m;
  l1c_int p = ax_funs.p;

  *NP = (l1c_NestaProb){.n=n, .m=m, .p=p, .fx=0.0, .xo=NULL, .xk=NULL, .yk=NULL, .zk=NULL,
                        .Atb=NULL, .gradf=NULL, .gradf_sum=NULL,
                        .dwork1=NULL, .dwork2=NULL, .b=NULL, .ax_funs={0},
                        .sigma=0, .mu=0, .tol=0, .L=0};


  NP->xo = l1c_calloc_double(m);
  NP->xk = l1c_calloc_double(m);
  NP->yk = l1c_calloc_double(m);
  NP->zk = l1c_calloc_double(m);
  NP->Atb = l1c_calloc_double(m);
  NP->gradf = l1c_calloc_double(m);
  NP->gradf_sum = l1c_calloc_double(m);
  NP->dwork1 = l1c_calloc_double(imax(m, p));
  NP->dwork2 = l1c_calloc_double(imax(m, p));

  if (!NP->xo ||!NP->xk || !NP->yk || !NP->zk
      || !NP->Atb || !NP->gradf || !NP->gradf_sum
      || !NP->dwork1 || !NP->dwork2){
    l1c_free_nesta_problem(NP);
    return NULL;
  }

  return NP;
}

/*
   Release the memory allocated l1c_init_nesta_problem().
   The memory for .b and .ax_funs is not released, since those were not allocated
   by l1c_init_nesta_problem().
 */
void l1c_free_nesta_problem(l1c_NestaProb *NP){
  if (NP){
    l1c_free_double(NP->xo);
    l1c_free_double(NP->xk);
    l1c_free_double(NP->yk);
    l1c_free_double(NP->zk);
    l1c_free_double(NP->Atb);
    l1c_free_double(NP->gradf);
    l1c_free_double(NP->gradf_sum);
    l1c_free_double(NP->dwork1);
    l1c_free_double(NP->dwork2);

    free(NP);
  }
}

/*
 * Implements the projections of yk and zk onto the (primal) feasible set Qp.
 * These are described in in eqs (3.5)-(3.7) and (3.10)-(3.12). Both can be put
 * the common framework of
 *
 * \f{align}{
 *  vk = \textrm{arg } \min_{x\in Q_p} \frac{L_{\mu}}{2} ||x - xx||_2^2 + \langle g, x\rangle
 * \f}
 *
 * Note that we have droped \f$x_k\f$ from the inner product, because it is a constant,
 * which justifies pushing (in step 3) the sum into the inner product, as is done in
 * in (3.10)
 *
 *
 * @note{By default, will use ax_funs.Ex() and ax_funs.Ety(), i.e., assuming the
 * analysis formulation. To use the synthesis formulation, set ax_funs.Ex=NULL
 * and ax_funs.Ety=NULL, and ax_funs.Ax will be used.}
 *
 * @param [in] NP
 * @param [in] xx
 * @param [in] g
 * @param [out] vk Solution vector
 */
void l1c_nesta_project(l1c_NestaProb *NP, double *xx, double *g, double *vk){

  l1c_int n = NP->n;
  l1c_int m = NP->m;

  double *Aq = NP->dwork2;
  double *AtAq = NP->dwork1;
  double Lmu = NP->L / NP->mu_j;
  double a0=0, a1=0, lambda=0;

  /* Store q in vk.*/
  l1c_daxpy_z(m, (-1.0/Lmu), g, xx, vk);

  NP->ax_funs.Ax(vk, Aq);
  NP->ax_funs.Aty(Aq, AtAq);

  // a0 = Lmu * [ ||Aq - b||/sigma  - 1]
  a0 = l1c_dnrm2_err(n, NP->b, Aq);
  a0 = Lmu * (a0/NP->sigma - 1.0);

  lambda = max(0, a0);

  a1 = lambda / (Lmu + lambda);


  /* We start with vk = q, and will compute
    vk = lambda/Lmu*(1 - a1)*Atb + q - a1*AtAq */

  cblas_daxpy(m, (lambda/Lmu) * (1.0 - a1), NP->Atb, 1, vk, 1);
  cblas_daxpy(m, -a1, AtAq, 1, vk, 1);

}


/*
 * Evaluate the (smoothed) functional and compute the gradient.
 *
 * @param [in,out] NP On exist, NP->gradf and NP->fx will be updated.

 y = Ax = Rx,
 y in R^n
 x in R^m
 u = W^T * x in R^p (as in ||W^T * x||_1)
 R (or A) is n by m
 W (or U) is p by m, u is in R^p

 See also eq (6.2). The transpose convention is different from their matlab code.
 */
void l1c_nesta_feval(l1c_NestaProb *NP){
  // l1c_int m = NP->m;
  l1c_int p = NP->p;
  double nrm_u2=0;
  double *Wtxk = NP->dwork1;
  double *u = NP->dwork2;

  /* If E and Et are void, we are doing synthesis, otherwise, analysis.
   */
  NP->ax_funs.Wtx(NP->xk, Wtxk);

  for (int i=0; i<p; i++){
    u[i] = Wtxk[i] / max(NP->mu_j, fabs(Wtxk[i]));
  }

  nrm_u2 = cblas_dnrm2(p, u, 1);
  nrm_u2 *= nrm_u2;

  NP->fx = cblas_ddot(p, u, 1, Wtxk, 1) - 0.5 * NP->mu_j * nrm_u2;

  NP->ax_funs.Wz(u, NP->gradf);
}


/*
 * Populates an l1c_NestaProb instance NP. NP should already be allocated by
 *  l1c_init_nesta_problem().
 *
 * @param [out] NP problem instance.
 * @param [in,out] beta_mu Initial beta such that `beta^n_con * mu0 = mu_final`
 * @param [in,out] beta_tol Factor such such that `beta_tol^n_con * tol0 = tol_final`
 * @param [in,out] b
 * @param [in] ax_funs
 * @param [in,out] opts
 *
 */
int l1c_nesta_setup(l1c_NestaProb *NP, double *beta_mu, double *beta_tol,
                    double *b, l1c_AxFuns ax_funs, l1c_NestaOpts *opts){

  double L=0;
  if (!ax_funs.Wz || !ax_funs.Wtx || !ax_funs.Rx || !ax_funs.Rty) {
    return L1C_INCONSISTENT_ARGUMENTS;
  }
  L = ax_funs.norm_W * ax_funs.norm_W;

  NP->b = b;
  NP->n_continue = opts->n_continue;
  NP->ax_funs = ax_funs;
  NP->sigma = opts->sigma;
  NP->mu = opts->mu;
  NP->L = L;
  NP->tol = opts->tol;

  double mu_final = opts->mu;
  double tol_final = opts->tol;
  double *Wtx_ref = NP->dwork1;

  /*Section 3.6, paragraph preceeding (3.14) suggests that
   \mu_0 = 0.9||A^Tb||_{\infty}, and we have A^T = W^T * R^T *b
  */
  NP->ax_funs.Aty(NP->b, NP->Atb);
  cblas_dcopy(NP->m, NP->Atb, 1, NP->xo, 1);

  NP->ax_funs.Wtx(NP->xo, Wtx_ref);
  l1c_abs_vec(NP->m, Wtx_ref, Wtx_ref);
  double mu0 = 0.9 * l1c_max_vec(NP->p, Wtx_ref);

  /** @todo tol=0.1 is a hueristic. What happens if user specifies tol > 0.1 ?
   */
  double tol0 = 0.1;

  /* We need mu_final and beta such that beta^n_cont * mu0= mu_final
     i.e., (beta^n_cont) = (mu_final/mu0)
     i.e., log(beta) = log(mu_final/mu0) / n_cont
     ie.,  beta = exp(log(mu_final/mu0) / n_cont)
   */

  *beta_mu = exp(log((mu_final / mu0)) / (double)opts->n_continue);
  *beta_tol = exp(log((tol_final / tol0)) / (double)opts->n_continue);

  /* After n continuation steps, NP->mu will again be what we said.*/
  NP->mu_j = mu0;
  NP->tol_j = tol0;

  return L1C_SUCCESS;

}

/** @ingroup nesta
 * Solves the problem
 *
 *   \f{align}{
 *   \min_x || W^T x||_1 \quad \textrm{s.t.}\quad ||Rx - b||_2 < \sigma
 *   \f}
 *
 *   if `opts.flags | L1C_ANALYSIS`
 *
 *   or
 *
 *   \f{align}{
 *   \min_x || x||_1 \quad \textrm{s.t.}\quad ||Ax - b||_2 < \sigma
 *   \f}
 *
 *   if `opts.flags | L1C_SYNTHESIS`
 *
 *   The matrix operators should be defined in ax_funs such that
 *
 *   `ax_funs.Mx =` \f$W\f$,  `ax_funs.Mty` = \f$ W^T\f$
 *   `ax_funs.Ex =` \f$Rx\f$, `ax_funs.Ety`= \f$ R^Ty\f$
 *
 *   FOr synthesis, `ax_funs.Mx`, and `ax_funs.Ex` may be `null` and
 *   we only need
 *   `ax_funs.Ax =` \f$Ax\f$, `ax_funs.Aty=` \f$A^Ty\f$
 *
 *   In general, `A` (or `R`) is `n` by `m`, `n < m` and `W` is `m` by `p`, with `p >=m`.
 *
 *   @param[in] m Length of true signal and xk.
 *   @param[out] xk On exit, contains optimal value.
 *   @param[in] n Length of the measurements b.
 *   @param[in] b Measurement vector.
 *   @param[in] ax_funs See above.
 *   @param[in] opts See l1c_NestaOpts.
 */
int l1c_nesta(l1c_int m, double *xk, l1c_int n, double *b,
              l1c_AxFuns ax_funs, l1c_NestaOpts opts){

  int status=0; //, idx_fmu=0;


  double alpha_k=0, tau_k = 0;
  double fbar=0, rel_delta_fmu;
  double beta_mu, beta_tol;


  l1c_NestaProb *NP = _l1c_NestaProb_new(ax_funs);

  if (!NP ){
    return L1C_OUT_OF_MEMORY;
  }

  /* Initialize*/
  status += l1c_nesta_setup(NP,  &beta_mu, &beta_tol, b, ax_funs, &opts);
  if (status){
    return status;
  }

  struct l1c_fmean_fifo fbar_fifo = _l1c_new_fmean_fifo();

  cblas_dcopy(NP->m, NP->xo, 1, NP->xk, 1);

  for (int iter=1; iter<= opts.n_continue; iter++){
    NP->mu_j = NP->mu_j * beta_mu;
    NP->tol_j = NP->tol_j * beta_tol;

    /* Reset everthing.*/
    l1c_init_vec(L1C_NESTA_NMEAN, fbar_fifo.f_vals, 0);
    fbar_fifo.next = fbar_fifo.f_vals;
    fbar_fifo.n_total = 0;

    l1c_init_vec(m,  NP->gradf_sum, 0.0);

    if (opts.verbose > 0) {
      l1c_printf("Starting nesta continuation iter %d, with muj = %f\n", iter,
             NP->mu_j);
      l1c_printf("Iter |     fmu     |  Rel. Vartn fmu |\n");
      l1c_printf("------------------------------------------------------\n");
    }
    /* ---------------------------- MAIN ITERATION -------------------------- */
    for (int k=0; k < L1C_NESTA_MAX_INNER_ITER; k++){
        l1c_nesta_feval(NP);

      /* ----------------- Update yk ------------------------- */
      l1c_nesta_project(NP, NP->xk, NP->gradf, NP->yk);

      /* From paragraph right before (2.3). Reference [43,50] proves convergence
         for these alpha_k, tau_k
      */
      alpha_k = 0.5 * (double)(k + 1);
      tau_k = 2.0 / (double)(k + 3);

      /* Update cummulative, weighted sum of gradients (3.8)-(3.12).*/
      cblas_daxpy(m, alpha_k, NP->gradf, 1, NP->gradf_sum, 1);

      /* Projection for zk */
      l1c_nesta_project(NP, NP->xo, NP->gradf_sum, NP->zk);

      /* ----------------- Update xk ----------------
         xk = tau_k * zk + (1-tau_k) * yk  */
      l1c_daxpby_z(m, tau_k, NP->zk, (1-tau_k), NP->yk, NP->xk);

      /*------------ Check for exit -----------------
        Must compute fbar first, because it should not include current fx.
       */

      fbar = _l1c_mean_fmean_fifo(&fbar_fifo);
      _l1c_push_fmeans_fifo(&fbar_fifo, NP->fx);
      rel_delta_fmu = fabs(NP->fx - fbar) / fbar;
      if (opts.verbose > 0 && (((k+1) % opts.verbose) == 0)) {
        l1c_printf("%d     %.3e       %.2e   \n", k+1, NP->fx, rel_delta_fmu);
      }

      if (rel_delta_fmu < NP->tol_j){
        if (opts.verbose > 0){
          l1c_printf("   stopping: delta_fmu < tol (%g < tol %g )\n\n", rel_delta_fmu, NP->tol_j);
        }
        break;
      }
    } /* Inner iter*/

    cblas_dcopy(NP->m, NP->xk, 1, NP->xo, 1);
  }

  cblas_dcopy(m, NP->xk, 1, xk, 1);
  l1c_free_nesta_problem(NP);
  free(fbar_fifo.f_vals);
  return status;
}
