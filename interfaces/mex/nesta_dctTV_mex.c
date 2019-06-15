#include "config.h"
/* It is important that mex.h is included before l1c_common.h
 so that, for this function, we get matlabs definition of fprintf
*/
#include "cblas.h"
#include "matrix.h"
#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "l1c.h"
#include "nesta.h"
#include "l1qc_newton.h"
#include "l1c_mex_utils.h"
#include "interfaces.h"


/*

 */
void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  l1c_NestaOpts opts = {.mu=1e-5,
                        .sigma=1e-3,
                        .tol=1e-3,
                        .n_continue=5,
                        .flags=L1C_SYNTHESIS};

  double *x_out=NULL,  *b=NULL;
  double *pix_idx_double=NULL, *x_ours=NULL;
  double alp_v = 0, alp_h = 0;

  double *dpar;
  int size_dpar;

  l1c_int n=0, mrow=0, mcol=1, mtot=0, idx=0;
  l1c_int size_pix_idx=0;

  l1c_int *pix_idx;
  int status=0;

  /* --------- Parse and verify the Input ------------ */
  _mex_assert_num_outputs(nlhs, 2);
  _mex_assert_num_inputs(nrhs, 5);

  mrow = (l1c_int)_mex_get_double_scalar_or_fail(prhs, 0);
  mcol = (l1c_int)_mex_get_double_scalar_or_fail(prhs, 1);
  mtot = mcol * mrow;

  /* -------- Check b -------------*/
  _mex_get_double_array_or_fail(prhs, 2, &b, &n);

  /* -------- Check pix_idx -------------*/
  _mex_get_double_array_or_fail(prhs, 3, &pix_idx_double, &size_pix_idx);

  /* -------- Check dpar -------------*/
  _mex_get_double_array_or_fail(prhs, 4, &dpar, &size_dpar);

  opts.verbose = dpar[NESTA_VERBOSE_IDX];
  opts.sigma = dpar[NESTA_SIGMA_IDX];
  opts.mu = dpar[NESTA_MU_IDX];
  opts.tol = dpar[NESTA_TOL_IDX];
  opts.n_continue = (int)dpar[NESTA_NCONT_IDX];
  opts.flags = (unsigned)dpar[NESTA_MODE_IDX];
  alp_v = dpar[NESTA_ALPV_IDX];
  alp_h = dpar[NESTA_ALPH_IDX];

  /* I think the problem is that we get the norm of U wrong when we try this,
     because in synthesis mode, ||U|| = 1.
     Probably can be fixed by passing flags to ax_funs, and adjusting how we
     set everything up. Then nesta can be dumb to which mode it is in as well.
   */
  if ( (alp_v >0 || alp_h > 0) && (opts.flags & L1C_SYNTHESIS)){
    mexErrMsgIdAndTxt("l1c:InvalidInput",
                      "Can only do TV with analysis mode.");
  }
  /* Ensure the sizes are consistent. */
  if ( !(mrow*mcol > n) || (n != size_pix_idx) || mrow <=0 || mcol <=0){
    mexErrMsgIdAndTxt("l1c:incompatible_dimensions",
                      "Must have length(x0) > length(b), and length(b) = length(pix_idx).");
  }
  if (opts.verbose > 0){
    printf("Input Parameters\n---------------\n");
    printf("   sigma:         %.5e\n", opts.sigma);
    printf("   mu:              %.5e\n", opts.mu);
    printf("   n_continue:      %d\n", opts.n_continue);
    printf("   alp_v:             %.5e\n", alp_v);
    printf("   alp_h:             %.5e\n", alp_h);
  }




  /* We are going to change x, so we must allocate and make a copy, so we
     dont change data in Matlabs workspace.
  */
  x_ours = l1c_malloc_double(mrow*mcol);
  pix_idx = calloc(n, sizeof(l1c_int));
  if (!x_ours || !pix_idx){
    status = L1C_OUT_OF_MEMORY;
    goto exit1;
  }
  /*
    pix_idx will naturally be a double we supplied from matlab
    and will have 1-based indexing. Convert to integers and
    shift to 0-based indexing and validate the indeces.  */
  for (int i=0; i<n; i++){
    idx = ((l1c_int) pix_idx_double[i]) - 1;
    if (idx <0 || idx > mtot-1){
      status = L1C_INVALID_ARGUMENT;
      goto exit1;
    }
    pix_idx[i] = idx;
  }
  l1c_AxFuns ax_funs;
  status |= l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, pix_idx, &ax_funs);
  if (status){
    goto exit1;
  }

  status |= l1c_nesta(mtot, x_ours, n, b, ax_funs, opts);


  /* Prepare output data.
   */
  plhs[0] = mxCreateDoubleMatrix((mwSize)mrow, (mwSize)mcol, mxREAL);
  plhs[1] = mxCreateDoubleScalar((double)status);

  x_out = mxGetPr(plhs[0]);
  cblas_dcopy(mrow * mcol, x_ours, 1, x_out, 1);
  // double *stat_out = mxGetPr(plhs[1]);
  // *stat_out = (double)status;

  ax_funs.destroy();

 exit1:
  free(pix_idx);
  l1c_free_double(x_ours);

  switch (status){
  case L1C_OUT_OF_MEMORY:
    mexErrMsgIdAndTxt("l1c:l1qc_dct:outofmemory",
                      "Error Allocating memory.");
    break;
  case L1C_INVALID_ARGUMENT:
    mexErrMsgIdAndTxt("l1c:l1qc_dct:bad_index",
                      "pix_idx contained invalid value(s).");
    break;

  case L1C_DCT_INIT_FAILURE:
    mexErrMsgIdAndTxt("l1c:dct_failure",
                      "Failed to initialize DCT transforms\n");
    break;
  default:
    break;
  }

} /* ------- mexFunction ends here ----- */
