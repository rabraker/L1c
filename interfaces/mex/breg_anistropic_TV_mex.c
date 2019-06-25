#include "config.h"
/* It is important that mex.h is included before l1c_common.h
 so that, for this function, we get matlabs definition of fprintf
*/
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "cblas.h"
#include "l1c.h"
#include "l1c_mex_utils.h"

#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

/**

 The matlab protype is
 [f_opt] = l1qc_dct(f, mu, tol, max_iter);
 where
 f an n by m matrix, n>2, m>2, mu and tol are scarlar doubles.
 */
void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  // breg_anistropix_TV(f, mu, tol, max_iter);

  double *f=NULL, *f_ours=NULL, *f_out=NULL, *uk=NULL;

  double tol=0, mu=0;
  l1c_int i=0, N=0, M=0;
  /*--- Defaults --- */
  int max_iter = 1000, max_jac_iter = 1;

  _mex_assert_num_outputs(nlhs, 1);
  _mex_assert_num_inputs(nrhs, 4);

  /* -------- Check f -------------*/
  _mex_assert_double(prhs, 0);
  /* check that f is a matrix, at least 3 by 3. */
  _mex_assert_2Darray_with_size(prhs, 0, 3, 3);

  /*Matlab give us column major order.*/
  M = (l1c_int) mxGetN(prhs[0]);
  N = (l1c_int) mxGetM(prhs[0]);
  f = mxGetPr(prhs[0]);

  /* -------- Check mu, tol, max_iter -------------*/
  mu = _mex_get_double_scalar_or_fail(prhs, 1);
  tol = _mex_get_double_scalar_or_fail(prhs, 2);
  max_iter = (int)_mex_get_double_scalar_or_fail(prhs, 3);

  /* We require that f is aligned on a DALIGN byte boundary. Matlab does not guarantee this.
  */
  f_ours = l1c_malloc_double(N*M);
  uk = l1c_malloc_double(N*M);
  if (!f_ours || !uk){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:outofmemory",
                      "Error Allocating memory.");
    goto exit;
  }
  cblas_dcopy(N*M, f, 1, f_ours, 1);
  int stat = l1c_breg_anistropic_TV(N, M, uk, f_ours, mu, tol, max_iter, max_jac_iter);
  if (stat){
    plhs[0] = NULL;
    goto exit;
  }
  /* Prepare output data. */
  plhs[0] = mxCreateDoubleMatrix((mwSize)M,  (mwSize)N, mxREAL);

  f_out =  mxGetPr(plhs[0]);

  for (i=0; i<N*M; i++){
    f_out[i] = uk[i];
  }

 exit:
  l1c_free_double(f_ours);
  l1c_free_double(uk);

} /* ------- mexFunction ends here ----- */
