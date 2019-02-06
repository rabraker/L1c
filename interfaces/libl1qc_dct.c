#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "l1c_common.h"
#include "l1qc_newton.h"

/* dct_mkl.h defines the Ax and Aty operations.
   To adapt this mex file to a different set of transformations,
   this largely what must be changed. Your new set of transformations
   must expose three functions with the following prototype:

   void Ax(double *x, double *y)
   void Aty(double *y, double *x)
   void AtAx(double *x, double *z)

   where x and z have length n, and y has length m and m<n.
*/

#include "dct_mkl.h"


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
/*
  [("epsilon", c_double),
  ("mu", c_double),
  ("lbtol", c_double),
  ("newton_tol", c_double),
  ("newton_max_iter", c_int),
  ("verbose", c_int),
  ("l1_tol", c_double)
  ("cgtol", c_double),
  ("cgmaxiter", c_int),
  ("warm_start_cg", c_int),
]
 */
int l1qc_dct(int N, double *x_out, int M, double *b, l1c_int *pix_idx,
             L1qcDctOpts opts, LBResult *lb_res){
  /* Ensure intel doesnt fuck us.*/
  mkl_set_interface_layer(MKL_INTERFACE_ILP64);
  mkl_set_threading_layer(MKL_THREADING_GNU);


  /*
    Pointers to the transform functions that define A*x and A^T *y
   */
  AxFuns Ax_funs = {.Ax=dctmkl_EMx_new,
                    .Aty=dctmkl_MtEty,
                    .AtAx=dctmkl_MtEt_EMx_new};
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
                       .warm_start_cg=opts.warm_start_cg};

  // LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};


  /* Allocate memory for x and b, which is aligned to DALIGN.
     Pointer address from caller probably wont be properly aligned.
   */
  printf("N = %d, M=%d\n", N, M);
  double *eta_0 = malloc_double(N);
  double *b_ours = malloc_double(M);

  if ( !b_ours || !eta_0){
    fprintf(stderr, "Memory Allocation failure\n");
    free_double(b_ours);
    free_double(eta_0);
    return 1;
  }

  cblas_dcopy(M, b, 1, b_ours, 1);

  dctmkl_setup(N, M, (l1c_int*)pix_idx);
  Ax_funs.Aty(b, eta_0);
  *lb_res = l1qc_newton(N, eta_0, M, b, params, Ax_funs);

  /* We solved for eta in the DCT domain. Transform back to
     standard coorbbdinates.
  */
  dctmkl_idct(eta_0);
  cblas_dcopy(N, eta_0, 1, x_out, 1);

  /* Cleanup our mess. */
  dctmkl_destroy();

  free_double(b_ours);
  free_double(eta_0);
  mkl_free_buffers();

  return 0;

}

