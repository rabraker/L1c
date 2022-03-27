#include "config.h"
#include <math.h>

#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#ifdef L1C_MEX_MATLAB
#include "matrix.h"
#endif
#include "mex.h"

#include "l1c.h"
#include "nesta.h"
#include "l1qc_newton.h"
#include "l1c_mex_utils.h"
#include "l1c_logging.h"

/*

 */
void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] ){
  /* Replace printf with mexPrintf */
  l1c_replace_printf(mexPrintf);

  l1c_NestaOpts opts = {.mu = 1e-5,
                        .sigma = 1e-3,
                        .tol = 1e-3,
                        .n_continue = 5};

  double *z_out = NULL, *b = NULL, * b_ours = NULL;
  double *pix_idx_double=NULL, *x0=NULL;
  double alp_v = 0, alp_h = 0;

  l1c_int n=0, mrow=0, mcol=1, mtot=0, idx=0;
  l1c_int size_pix_idx=0;
  DctMode dct_mode;
  BpMode  bp_mode;
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
  _mex_assert_scalar_struct(prhs, 4);

  opts.sigma      = _mex_get_double_from_struct_or_fail(prhs, 4, "sigma");
  opts.verbose    = _mex_get_double_from_struct_or_fail(prhs, 4, "verbose");
  opts.mu         = _mex_get_double_from_struct_or_fail(prhs, 4, "mu");
  opts.tol        = _mex_get_double_from_struct_or_fail(prhs, 4, "tol");
  opts.n_continue = (int)_mex_get_double_from_struct_or_fail(prhs, 4, "n_continue");
  alp_v           = _mex_get_double_from_struct_or_fail(prhs, 4, "alpha_v");
  alp_h           = _mex_get_double_from_struct_or_fail(prhs, 4, "alpha_h");
  dct_mode        = (DctMode)_mex_get_double_from_struct_or_fail(prhs, 4, "dct_mode");
  mxArray *mxA_mode = mxGetField(prhs[4], 0, "bp_mode");
  if (!mxA_mode) {
    mexErrMsgIdAndTxt("l1c:notAField", "Error loading field '%s'.", "mode");
  }
  // [synthesis | analysis]
  char mode[10];
  if (mxGetString(mxA_mode, mode, 10)){
    mexErrMsgIdAndTxt("l1c:badField", "Field 'mode' of opts must be a string, either analysis or synthesis");
  }

  if (strncmp(mode, "synthesis", 10) == 0){
    bp_mode = synthesis;
  }else if(strncmp(mode, "analysis", 10) == 0){
    bp_mode = analysis;
  }else{
    mexErrMsgIdAndTxt("l1c:badField", "Field 'mode' of opts must be a string, either analysis or synthesis");
    /*We never reach this return, but it prevents warnings about unitialized bp_mode.*/
    return;
  }


  /* Ensure the sizes are consistent. */
  if ( !(mrow*mcol > n) || (n != size_pix_idx) || mrow <=0 || mcol <=0){
    mexErrMsgIdAndTxt("l1c:incompatible_dimensions",
                      "Must have length(x0) > length(b), and length(b) = length(pix_idx).");
  }
  if (opts.verbose > 0){
    l1c_printf("Input Parameters\n---------------\n");
    l1c_printf("   sigma:         %.5e\n", opts.sigma);
    l1c_printf("   mu:              %.5e\n", opts.mu);
    l1c_printf("   n_continue:      %d\n", opts.n_continue);
    l1c_printf("   alp_v:             %.5e\n", alp_v);
    l1c_printf("   alp_h:             %.5e\n", alp_h);
  }


  /* We are going to change x, so we must allocate and make a copy, so we
     dont change data in Matlabs workspace.
  */
  pix_idx = calloc(n, sizeof(l1c_int));

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
  status |= l1c_setup_dctTV_transforms(n, mrow, mcol, alp_v, alp_h, dct_mode,
                                       bp_mode, pix_idx, &ax_funs);
  if (status){
    goto exit1;
  }

  x0 = l1c_malloc_double(ax_funs.m);
  b_ours = l1c_malloc_double(n);
  if (NULL == pix_idx || NULL == b_ours || NULL == x0) {
    status = L1C_OUT_OF_MEMORY;
    goto exit1;
  }

  cblas_dcopy(n, b, 1, b_ours, 1);

  status |= l1c_nesta(mtot, x0, n, b_ours, ax_funs, opts);


  /* Prepare output data.
   */
  plhs[0] = mxCreateDoubleMatrix((mwSize)mrow, (mwSize)mcol, mxREAL);
  plhs[1] = mxCreateDoubleScalar((double)status);

  z_out = mxGetPr(plhs[0]);

  ax_funs.Mx(x0, z_out); // A straigt copy in analysis mode


  ax_funs.destroy();

 exit1:
  free(pix_idx);
  l1c_free_double(x0);
  l1c_free_double(b_ours);

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
