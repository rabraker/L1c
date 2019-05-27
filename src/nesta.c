#include <stdlib.h>
#include<cblas.h>
#include <math.h>
#include "l1c.h"
#include "l1c_math.h"
#include "nesta.h"
#define N_MEAN 10




/**
   Initialize a l1c_NestaProblem struct. Memory will be allocated for all
   arrays.
 */
l1c_NestaProb* l1c_init_nesta_problem(l1c_int n, l1c_int m, double *b, l1c_AxFuns ax_funs,
                                      double sigma, double mu, double tol, double L, unsigned flags){

  l1c_NestaProb *NP = malloc(sizeof(l1c_NestaProb));
  if (!NP){
    return NULL;
  }

  *NP = (l1c_NestaProb){.n=n, .m=m, .fx=0.0, .xo=NULL, .xk=NULL, .yk=NULL, .zk=NULL,
                    .Atb=NULL, .gradf=NULL, .gradf_sum=NULL,
                    .dwork1=NULL, .dwork2=NULL, .b=b, .ax_funs=ax_funs,
                        .sigma=sigma, .mu=mu, .tol=tol, .L_by_mu=L/mu, .flags=flags};


  NP->xo = l1c_calloc_double(m);
  NP->xk = l1c_calloc_double(m);
  NP->yk = l1c_calloc_double(m);
  NP->zk = l1c_calloc_double(m);
  NP->Atb = l1c_calloc_double(m);
  NP->gradf = l1c_calloc_double(m);
  NP->gradf_sum = l1c_calloc_double(m);
  NP->dwork1 = l1c_calloc_double(m);
  NP->dwork2 = l1c_calloc_double(m);

  if (!NP->xo ||!NP->xk || !NP->yk || !NP->zk
      || !NP->Atb || !NP->gradf || !NP->gradf_sum
      || !NP->dwork1 || !NP->dwork2){
    l1c_free_nesta_problem(NP);
    return NULL;
  }

  if(flags | L1C_ANALYSIS){
    NP->ax_funs.Ety(b, NP->Atb);
  }else{
    ax_funs.Aty(b, NP->Atb);
  }


  return NP;
}

/**
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

/**
 * By default, will use ax_funs.Ex() and ax_funs.Ety(), i.e., assuming the
 * analysis formulation. To use the synthesis formulation, set ax_funs.Ex=NULL
 * and ax_funs.Ety=NULL.
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
  double Lmu = NP->L_by_mu;
  double a0=0, a1=0, lambda=0;

  /* Store q in vk.*/
  l1c_daxpy_z(m, (-1.0/Lmu), g, xx, vk);

  if(NP->flags | L1C_ANALYSIS){
    NP->ax_funs.Ex(vk, Aq);
    NP->ax_funs.Ety(Aq, AtAq);
  }else{
    NP->ax_funs.Ax(vk, Aq);
    NP->ax_funs.Aty(Aq, AtAq);
  }

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


/**
 * Evaluate the (smoothed) functional and compute the gradient.
 *
 * @param [in,out] NP On exist, NP->gradf and NP->fx will be updated.
 */
void l1c_nesta_feval(l1c_NestaProb *NP){
  l1c_int m = NP->m;

  double nrm_u2=0;
  double *Uxk = NP->dwork1;
  double *u = NP->dwork2;

  /* If E and Et are void, we are doing synthesis, otherwise, analysis.
   */
  if( NP->flags | L1C_ANALYSIS){
    NP->ax_funs.Mx(NP->xk, Uxk);
  }else{
    cblas_dcopy(m, NP->xk, 1, Uxk, 1);
  }

  for (int i=0; i<m; i++){
    u[i] = Uxk[i] / max(NP->mu, fabs(Uxk[i]));
  }


  nrm_u2 = cblas_dnrm2(m, u, 1);
  nrm_u2 *= nrm_u2;

  NP->fx = cblas_ddot(m, u, 1, Uxk, 1) - 0.5 * NP->mu * nrm_u2;

  NP->ax_funs.Mty(u, NP->gradf);

}



int l1c_nesta(l1c_int m, double *xk, double mu, l1c_int n, double *b,
              l1c_AxFuns ax_funs, double sigma){

  int max_iter = 10000;
  int status=0; //, idx_fmu=0;


  double alpha_k=0, tau_k = 0;

  /*Lipschitz constant divided by mu */
  double L = 1;
  double tol = 1e-4;

  unsigned flags = L1C_ANALYSIS;

  // double fmu_store[N_MEAN+1];
  // for (int i=0; i<N_MEAN+1; i++){
  //   fmu_store[i] = 0;
  // }

  l1c_NestaProb *NP = l1c_init_nesta_problem(n, m, b, ax_funs, sigma, mu, tol, L, flags);

  if (!NP ){
    return L1C_OUT_OF_MEMORY;
  }


  l1c_init_vec(m,  NP->gradf_sum, 0.0);

  /* ---------------------------- MAIN ITERATION -------------------------- */
  for (int k=0; k < max_iter; k++){

    l1c_nesta_feval(NP);

    /* ----------------- Update yk ------------------------- */
    l1c_nesta_project(NP, NP->xk, NP->gradf, NP->yk);

    /* From paragraph right before (2.3). Reference [43,50] proves convergence
     for these alpha_k, tau_k
    */
    alpha_k = 0.5 * (double)(k + 1);
    tau_k = 2.0 / ((double)k + 3.0);

    /* Update cummulative, weighted sum of gradients (3.8)-(3.12).*/
    cblas_daxpy(m, alpha_k, NP->gradf, 1, NP->gradf_sum, 1);

    /* Projection for zk */
    l1c_nesta_project(NP, NP->xo, NP->gradf_sum, NP->zk);


    /* ---- Update xk ----------------*/
    /* xk = tau_k * zk + (1-tau_k) * yk*/
    l1c_daxpby_z(m, tau_k, NP->zk, (1-tau_k), NP->yk, NP->xk);
  }

  cblas_dcopy(m, NP->xk, 1, xk, 1);
  return status;
}
