#include "config.h"
/* It is important that mex.h is included before l1c_common.h
 so that, for this function, we get matlabs definition of fprintf
*/
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include "l1c_common.h"
#include "l1c_memory.h"
#include <math.h>
#include "l1c.h"


#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

/*
 *	m e x F u n c t i o n
 int breg_anistropic_TV(l1c_int n, l1c_int m, double *uk, double *f,
 double mu, double tol, int max_iter, int max_jac_iter);

 The matlab protype is
 [f_opt] = l1qc_dct(f, mu, tol, max_iter, max_jac_iter);
 where
 b has length m
 pix_idx has length m
 and opts is an options struct. See the l1qc_dct.m doc string for details.
 */
void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  // l1qc(x0, b, pix_idx, params);
  /* inputs */


  double *f=NULL, *f_ours=NULL, *f_out=NULL, *uk=NULL;

  /*--- Defaults --- */
  int max_iter = 1000, max_jac_iter = 1;
  double tol=0.001, mu = 5;

  l1c_int i=0, N=0, M=0;
  // mwSize *dims;
  if( !(nlhs > 0)) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:nlhs",
                      "One output required.");
  }

  if( (nrhs != 4) ){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:nrhs",
                      "four inputs required.");
  }

  /* -------- Check f -------------*/
  if( !mxIsDouble(prhs[0]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "First Input vector must be type double.");
  }

  /* check that f is a matrix, at least 3 by 3. */
  if( (mxGetN(prhs[2]) > 2)  && (mxGetM(prhs[2]) > 2) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notVector",
                      "Second Input must be at least a 3 by 3 matrix .");
  }else{
    /*Check, but I believe matlab give us column major order.*/
    M = (l1c_int) mxGetN(prhs[0]);
    N = (l1c_int) mxGetM(prhs[0]);
  }

  f = mxGetPr(prhs[0]);


  /* -------- Check mu -------------*/
  if( !mxIsDouble(prhs[1])) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "mu a scalar double.");
  }else if( nrhs > 1){
    mu = (double)mxGetScalar(prhs[1]);
  }

  /* -------- Check tol -------------*/
  if( !mxIsDouble(prhs[2])) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "mu a scalar double.");
  }else if(nrhs > 2){
    tol = (double)mxGetScalar(prhs[2]);
  }

  /* -------- Check max_iter -------------*/
  if( !mxIsDouble(prhs[3]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "mu a scalar double.");
  }else if(nrhs > 3) {
    max_iter = (int)mxGetScalar(prhs[3]);
  }


  // printf("Input Parameters\n---------------\n");
  // printf("   mu:           %f\n", mu);
  // printf("   tol:          %.5e\n", tol);
  // printf("   max_iter:     %d\n", max_iter);
  // printf("   max_iter_jac: %d\n", max_jac_iter);


  /* We require that f is aligned on a DALIGN byte boundary. Matlab does not guarantee this.
  */
  f_ours = malloc_double(N*M);
  uk = malloc_double(N*M);
  if (!f_ours || !uk){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:outofmemory",
                      "Error Allocating memory.");
  }
  cblas_dcopy(N*M, f, 1, f_ours, 1);
  int stat = breg_anistropic_TV(N, M, uk, f_ours, mu, tol, max_iter, max_jac_iter);
  if (stat){
    plhs[0] = NULL;
    goto exit;
  }
  /* Prepare output data.
   */
  plhs[0] = mxCreateDoubleMatrix((mwSize)M,  (mwSize)N, mxREAL);

  f_out =  mxGetPr(plhs[0]);

  for (i=0; i<N*M; i++){
    f_out[i] = uk[i];
  }

 exit:
 free_double(f_ours);
 free_double(uk);
 if (stat){
   mexErrMsgIdAndTxt("l1c:l1qc_dct:outofmemory",
                     "Error Allocating memory.");
 }
#if defined(_USEMKL_)
 mkl_free_buffers();
#endif
} /* ------- mexFunction ends here ----- */
