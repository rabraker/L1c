#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "l1c_common.h"
#include "l1qc_newton.h"
#include "omp.h"
#include <sys/time.h>
#include <unistd.h>

/* dct_mkl.h defines the Ax and Aty operations.
   To adapt this mex file to a different set of transformations,
   this largely what must be changed. Your new set of transformations
   must expose three functions with the following prototype:

   void Ax(double *x, double *y)
   void Aty(double *y, double *x)
   void AtAx(double *x, double *z)

   where x and z have length n, and y has length m and m<n.
*/

#ifdef _USEMKL_
#include "dct_mkl.h"

#define Ax_fun dctmkl_EMx_new
#define Aty_fun dctmkl_MtEty
#define AtAx_fun dctmkl_MtEt_EMx_new
#define ax_idct dctmkl_idct
#define ax_setup dctmkl_setup
#define ax_destroy dctmkl_destroy

#else
#ifdef _USEFFTW3_

#include "dct.h"
#define Ax_fun dct_EMx_new
#define Aty_fun dct_MtEty
#define AtAx_fun dct_MtEt_EMx_new
#define ax_idct dct_idct
#define ax_setup dct_setup
#define ax_destroy dct_destroy

#endif
#endif



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




int l1qc_dct(int N, double *x_out, int M, double *b, l1c_int *pix_idx,
             L1qcDctOpts opts, LBResult *lb_res){
  int status = 0;

  long start, end;
  struct timeval timecheck;
  gettimeofday(&timecheck, NULL);
  start = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec/1000;

#if defined(HAVE_LIBMKL_RT)
  /* Ensure intel doesnt fuck us.*/
  mkl_set_interface_layer(MKL_INTERFACE_ILP64);
  mkl_set_threading_layer(MKL_THREADING_GNU);
#endif

  /*
    Pointers to the transform functions that define A*x and A^T *y
  */
  AxFuns Ax_funs = {.Ax=Ax_fun,
                    .Aty=Aty_fun,
                    .AtAx=AtAx_fun};

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
    status =  L1C_OUT_OF_MEMORY;
    goto fail;
  }

  cblas_dcopy(M, b, 1, b_ours, 1);

  if (ax_setup(N, M, (l1c_int*)pix_idx)){
    status = L1C_DCT_INIT_FAILURE;
    goto fail;
  }
  Ax_funs.Aty(b, eta_0);

  *lb_res = l1qc_newton(N, eta_0, M, b, params, Ax_funs);

  /* We solved for eta in the DCT domain. Transform back to
     standard coorbbdinates.
  */
  ax_idct(eta_0);

  cblas_dcopy(N, eta_0, 1, x_out, 1);




  gettimeofday(&timecheck, NULL);
  end = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec/1000;

  double time_total = ((double)(end -start)) / 1000.0;

  printf("total c time: %f\n", time_total);
  printf("total-newton-iter: %d\n", lb_res->total_newton_iter);
  printf("total-cg-iter: %d\n", lb_res->total_cg_iter);
  printf("time per cg iter: %g\n", time_total / (double) lb_res->total_cg_iter);


  ax_destroy(); // Should not call this if ax_setup() failed.

 fail:
  /* Cleanup our mess. */
  free_double(b_ours);
  free_double(eta_0);
#ifdef _USEMKL_
  mkl_free_buffers();
#endif

  return status;
}

