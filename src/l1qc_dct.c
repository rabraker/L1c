#include "config.h"
#include "l1c_timing.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "l1c_common.h"
#include "l1c_math.h"
#include "l1c_memory.h"
#include "l1c_transforms.h"
#include "l1qc_newton.h"
#include "omp.h"


typedef struct L1qcDctOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  int lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cgtol;
  int cgmaxiter;
  int warm_start_cg;
}L1qcDctOpts;




/*--------------------------- Two-dimensional version ------------------------------*/

int l1qc_dct(int Nrow, int Ncol, double *x_out, int M, double *b, l1c_int *pix_idx,
             L1qcDctOpts opts, LBResult *lb_res){
  struct timeval tv_start, tv_end;
  tv_start = l1c_get_time();

  int status = 0;
  int Ntot = Nrow*Ncol;
  /*
    Struct of pointers to the transform functions that define A*x and A^T *y
  */
  L1cAxFuns ax_funs;

  NewtParams params = {.epsilon = opts.epsilon,
                       .tau = opts.tau,
                       .mu = opts.mu,
                       .newton_tol = opts.newton_tol,
                       .newton_max_iter = opts.newton_max_iter,
                       .lbiter = 0,
                       .l1_tol = opts.l1_tol,
                       .lbtol = opts.lbtol,
                       .verbose = opts.verbose,
                       .cg_params.max_iter = opts.cgmaxiter,
                       .cg_params.tol = opts.cgtol,
                       .cg_params.verbose=0,
                       .warm_start_cg=opts.warm_start_cg};


  /* Allocate memory for x and b, which is aligned to DALIGN.
     Pointer address from caller probably wont be properly aligned.
   */
  printf("Nrow = %d, Ncol=%d, M=%d\n", Nrow, Ncol, M);

  status = l1c_setup_dct_transforms(Nrow, Ncol, M, pix_idx, &ax_funs);

  if (status != 0 || !ax_funs.Ax || !ax_funs.Aty || !ax_funs.AtAx){
    fprintf(stderr, "Error setup up dct. Exiting.\n");
    status = L1C_DCT_INIT_FAILURE;
    return status;
  }

  double *eta_0 = malloc_double(Ntot);
  double *b_ours = malloc_double(M);

  if ( !b_ours || !eta_0){
    fprintf(stderr, "Memory Allocation failure\n");
    status =  L1C_OUT_OF_MEMORY;
    goto exit1;
  }

  /* eta_0 init is only needed for when we use MKL.*/
  l1c_init_vec(Ntot, eta_0, 0);
  cblas_dcopy(M, b, 1, b_ours, 1);

  ax_funs.Aty(b_ours, eta_0);

  *lb_res = l1qc_newton(Ntot, eta_0, M, b_ours, params, ax_funs);
  if (lb_res->status){
    status = lb_res->status;
    goto exit2;
  }
  /* We solved for eta in the DCT domain. Transform back to
     standard coorbbdinates.
  */
  ax_funs.M(eta_0);

  cblas_dcopy(Ntot, eta_0, 1, x_out, 1);


  tv_end = l1c_get_time();

  double time_total = l1c_get_time_diff(tv_start, tv_end);

  printf("total c time: %f\n", time_total);
  printf("total-newton-iter: %d\n", lb_res->total_newton_iter);
  printf("total-cg-iter: %d\n", lb_res->total_cg_iter);
  printf("time per cg iter: %g\n", time_total / (double) lb_res->total_cg_iter);

 exit2:
  ax_funs.destroy(); // Should not call this if ax_setup() failed.
 exit1:
  /* Cleanup our mess. */
  free_double(b_ours);
  free_double(eta_0);
  return status;
}

