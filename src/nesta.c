#include<cblas.h>
#include <math.h>
#include "l1c.h"
#include "l1c_math.h"

#define N_MEAN 10


void nesta_project(l1c_int m, double *x, double *g, double *Atb, double *vk, double *dwork1,
                   l1c_int n, double *b, double *dwork2,
                   double L_by_mu, l1c_AxFuns ax_funs, double sigma){

  double *Aq = dwork2;
  double *AtAq = dwork1;

  /* Store q in vk.*/

  double tmp2=0, lambda=0, tmp1=0;

  l1c_daxpy_z(m, (-1/L_by_mu), g, x, vk);
  if( ax_funs.E  && ax_funs.Et){
    ax_funs.E(vk, Aq);
    ax_funs.Et(Aq, AtAq);
  }else{
    ax_funs.Ax(vk, Aq);
    ax_funs.Aty(Aq, AtAq);

  }

  tmp1 = L_by_mu * l1c_dnrm2_err(n, b, Aq);
  tmp1 = tmp1/sigma - 1.0;

  lambda = max(0, tmp1);

  tmp2 = lambda / (lambda + L_by_mu);
  // vk = lambda/L_by_mu*(1 - gamma)*Atb + q - gamma*AtAq;

  cblas_daxpy(m, tmp1, Atb, 1, vk, 1);

  cblas_daxpy(m, -tmp2, AtAq, 1, vk, 1);

}


void nesta_l1_feval(l1c_int m, double *x, double *u, double *gradf, double *fx,
                    double mu, l1c_AxFuns ax_funs){

  double nrm_u2=0;
  /* If E and Et are void, we are doing synthesis, otherwise, analysis.
   */
  if( !ax_funs.E  && !ax_funs.Et){
    cblas_dcopy(m, x, 1, u, 1);
  }else{
    cblas_dcopy(m, x, 1, u, 1);
    ax_funs.M(u);
  }

  for (int i=0; i<m; i++){
    u[i] = u[i] / max(mu, fabs(u[i]));
  }


  nrm_u2 = cblas_dnrm2(m, u, 1);
  nrm_u2 *= nrm_u2;

  *fx = cblas_ddot(m, u, 1, x, 1) - 0.5*mu * nrm_u2;

  cblas_dcopy(m, u, 1, gradf, 1);
  ax_funs.Mt(gradf);
}


int l1c_nesta(l1c_int m, double *xk, double *x0, double mu, l1c_int n, double *b,
              l1c_AxFuns ax_funs, double sigma){

  int max_iter = 10000;
  int status=0; //, idx_fmu=0;
  double *Atb=NULL, *AtAAtb=NULL, *dwork1=NULL, *dwork2=NULL, *gradf_sum=NULL;
  double *gradf=NULL, *uk=NULL, *yk=NULL, *zk=NULL;
  double fx=0;

  double alpha_k=0, tau_k = 0;

  /*Lipschitz constant divided by mu */
  double L_by_mu = 1; /* ----- FIX THIS ----- */

  // double fmu_store[N_MEAN+1];
  // for (int i=0; i<N_MEAN+1; i++){
  //   fmu_store[i] = 0;
  // }

  Atb = l1c_calloc_double(m);
  AtAAtb = l1c_calloc_double(m);
  dwork1 = l1c_calloc_double(m);
  dwork2 = l1c_calloc_double(m);
  gradf = l1c_calloc_double(m);
  gradf_sum = l1c_calloc_double(m);
  uk = l1c_calloc_double(m);
  yk = l1c_calloc_double(m);
  zk = l1c_calloc_double(m);

  if (!Atb || !AtAAtb || !dwork1 || !dwork2 || !gradf_sum
      || !gradf || !uk || !yk || !zk){
    status = L1C_OUT_OF_MEMORY;
    goto exit;
  }


  /* ---------------------------- MAIN ITERATION -------------------------- */
  for (int k=0; k < max_iter; k++){


    nesta_l1_feval(m, xk, uk, gradf, &fx, mu, ax_funs);

    /* ----------------- Update yk ------------------------- */
    nesta_project(m, xk, gradf, Atb, yk, dwork1, n, b, dwork2, L_by_mu, ax_funs, sigma);

    /* From paragraph right before (2.3). Reference [43,50] proves convergence
     for these alpha_k, tau_k
    */
    alpha_k = 0.5 * (double)(k + 1);
    tau_k = 2.0 / ((double)k + 3.0);

    /* Update cummulative, weighted sum of gradients (3.8)-(3.12).*/
    cblas_daxpy(m, alpha_k, gradf, 1, gradf_sum, 1);

    /* Projection for zk */
    nesta_project(m, x0, gradf_sum, Atb, zk, dwork1, n, b, dwork2, L_by_mu, ax_funs, sigma);


    /* ---- Update xk ----------------*/
    /* xk = tau_k * zk + (1-tau_k) * yk*/
    l1c_daxpby_z(m, tau_k, zk, (1-tau_k), yk, xk);
  }


 exit:
  l1c_free_double(Atb);
  l1c_free_double(AtAAtb);
  l1c_free_double(dwork1);
  l1c_free_double(dwork2);
  l1c_free_double(gradf_sum);
  l1c_free_double(gradf);
  l1c_free_double(uk);
  l1c_free_double(yk);
  l1c_free_double(zk);

  return status;
}
