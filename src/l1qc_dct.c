#include "config.h"
#include "l1c_timing.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "omp.h"

#include "cblas.h"
#include "l1c.h"
#include "l1c_math.h"
#include "l1qc_newton.h"


/**
 * Interface to l1qc_newton using the DCT1/DCT2 transforms \f$Ax=b\f$
 *
 * The size of `A` is `n` by `mrow*mcol`.
 *
 * @param[in] mrow number of rows of the underlying signal `x`
 * @param[in] mcol number of columns of the underlying signal `x`
 * @param[out] x_out An array of size mrow*mcol, which will contain the
 *             result of the optimization. Assumed aligned to DALIGN boundary.
 * @param[in] n The size of the observation array, b. In general, n << mrow*mcol.
 * @param[in] b The array of measurements or observations. Need NOT be aligned.
 * @param[in] pix_idx An array of indeces, such that, given the true mcol*mrow vectorx
 *            then x[pix_idx] = b, in python notation.
 * @param[in] opts A struct of options for the optimization.
 * @param[out] lb_res A struct containing the information about
 *             the results of the optimization.
*/
int l1qc_dct(int mrow, int mcol, double *x_out, int n, double *b, l1c_int *pix_idx,
             l1c_L1qcOpts opts, l1c_LBResult *lb_res){
  struct timeval tv_start, tv_end;
  tv_start = l1c_get_time();

  int status = 0;
  int mtot = mrow*mcol;
  /*
    Struct of pointers to the transform functions that define A*x and A^T *y
  */
  l1c_AxFuns ax_funs;


  /* Allocate memory for x and b, which is aligned to DALIGN.
     Pointer address from caller probably wont be properly aligned.
   */
  printf("mrow = %d, mcol=%d, n=%d\n", mrow, mcol, n);

  status = l1c_setup_dct_transforms(n, mrow, mcol, pix_idx, &ax_funs);

  if (status != 0 || !ax_funs.Ax || !ax_funs.Aty || !ax_funs.AtAx){
    fprintf(stderr, "Error setup up dct. Exiting.\n");
    status = L1C_DCT_INIT_FAILURE;
    return status;
  }

  double *eta_0 = l1c_malloc_double(mtot);
  double *b_ours = l1c_malloc_double(n);

  if ( !b_ours || !eta_0){
    fprintf(stderr, "Memory Allocation failure\n");
    status =  L1C_OUT_OF_MEMORY;
    goto exit1;
  }

  /* eta_0 init is only needed for when we use MKL.*/
  l1c_init_vec(mtot, eta_0, 0);
  cblas_dcopy(n, b, 1, b_ours, 1);

  ax_funs.Aty(b_ours, eta_0);

  *lb_res = l1c_l1qc_newton(mtot, eta_0, n, b_ours, opts, ax_funs);
  if (lb_res->status){
    status = lb_res->status;
    goto exit2;
  }
  /* We solved for eta in the DCT domain. Transform back to
     standard coorbbdinates.
  */

  ax_funs.Mx(eta_0, x_out);

  cblas_dcopy(mtot, eta_0, 1, x_out, 1);


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
  l1c_free_double(b_ours);
  l1c_free_double(eta_0);
  return status;
}

