#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "l1qc_newton.h"
#include "l1qc_common.h"
/*
 *	m e x F u n c t i o n
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  // l1qc(x0, b, pix_idx, params);
  /* inputs */
  double *x=NULL, *u=NULL, *b=NULL;
  NewtParams params = {.epsilon=0, .tau=0, .mu=0,
                        .newton_tol=0, .newton_max_iter = 0, .lbiter=0,
                        .lbtol=0, .verbose = 0, .cg_params.tol=0};
  int i=0,N=0, M=0, npix=0;
  double *pix_idx_=NULL, *x_out=NULL;
  int *pix_idx;
  // mwSize *dims;

  if(nrhs != 4) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:nrhs",
                      "four inputs required.");
  }

  if(nlhs != 1) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:nlhs",
                      "One output required.");
  }

  /* -------- Check x0 -------------*/
  if( !mxIsDouble(prhs[0]) ) {
      mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "First Input vector must be type double.");
  }

  /* check that number of rows in second input argument is 1 */
  // dims = mxGetDimensions(prhs[1]);
  if( (mxGetN(prhs[0]) > 1) ){
    N = mxGetN(prhs[0]);
  }else if( (mxGetM(prhs[0]) >1) ) {
    N = mxGetM(prhs[0]);
  }else{
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notRowVector",
                      "First Input must be a vector.");
  }

  /* -------- Check b -------------*/
  if( !mxIsDouble(prhs[1]) ) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "First Input vector must be type double.");
  }

  /* check that number of rows in second input argument is 1 */
 if( mxGetN(prhs[1]) > 1 ){
    M = mxGetN(prhs[1]);
  }else if( mxGetM(prhs[1]) >1 ) {
    M = mxGetM(prhs[1]);
  }else{
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notRowVector",
                      "Second Input must be a vector.");
  }

 //  /* -------- Check pix_idx -------------*/
  if( !mxIsDouble(prhs[2]) ) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "pix_idx vector must be type double.");
  }

  /* check that number of rows in second input argument is 1 */
 if( mxGetN(prhs[2]) > 1 ){
    npix = mxGetN(prhs[2]);
  }else if( mxGetM(prhs[2]) > 1 ) {
    npix = mxGetM(prhs[2]);
  }else{
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notRowVector",
                      "pix_idx must be a vector.");
  }

 if ( !(N > M) | (M != npix) ){
   printf("N = %d, M=%d, npix = %d\n", N, M, npix);
   mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notRowVector",
                     "Must have length(x0) > length(b), and length(b) = length(pix_idx).");
 }

 if ( !mxIsStruct(prhs[3])){
   mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notStruct",
                     "params must be a struct.");
 }
  x = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);


  pix_idx_ = mxGetPr(prhs[2]);
  pix_idx = calloc(M, sizeof(int));
  for (i=0; i<M; i++){
    pix_idx[i] = (int) pix_idx_[i];
  }

  printf("N = %d, M = %d, npix = %d\n", N, M, npix);


  int nfld = mxGetNumberOfFields(prhs[3]);
  // char *flds[]={ "verbose", "tau", "mu"};
  const char *name;

  mxArray *tmp;
  for (i=0; i<nfld; i++){
    // tmp = mxGetField(prhs[3], 0, "verbose");
    tmp = mxGetFieldByNumber(prhs[3], 0, i);
    name =mxGetFieldNameByNumber(prhs[3], i);
    if (!tmp){
      mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:norverbose",
                        "Error loading field '%s'.", name);
    }else if( !mxIsScalar(tmp) | !mxIsNumeric(tmp) ){
        mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:norverbose",
                          "Bad data in field '%s'. Fields in params struct must be numeric scalars.", name);
    }

    printf("name = %s, value=%f\n", name, mxGetScalar(tmp));
    if ( strcmp(name, "epsilon") == 0){
      params.epsilon = mxGetScalar(tmp);
    }else if ( strcmp(name, "mu") == 0){
      params.mu = mxGetScalar(tmp);
    }else if ( strcmp(name, "tau") == 0){
      params.tau = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_tol") == 0){
      params.newton_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_max_iter") == 0){
      params.newton_max_iter = (int) mxGetScalar(tmp);
    }else if ( strcmp(name, "lbiter") == 0){
      params.lbiter = (int)mxGetScalar(tmp);
    }else if ( strcmp(name, "lbtol") == 0){
      params.lbtol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgtol") == 0){
      params.cg_params.tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgmaxiter") == 0){
      params.cg_params.max_iter = mxGetScalar(tmp);
    }else if ( strcmp(name, "verbose") == 0){
      params.verbose= (int) mxGetScalar(tmp);
    }else{
      mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:norverbose",
                        "Unrecognized field '%s' in params struct.", name);
    }

  }

  printf("Input Parameters\n---------------\n");
  printf("   verbose:         %d\n", params.verbose);
  printf("   epsilon:         %.5e\n", params.epsilon);
  printf("   tau:             %.5e\n", params.tau);
  printf("   mu:              %.5e\n", params.mu);
  printf("   newton_tol:      %.5e\n", params.newton_tol);
  printf("   newton_max_iter: %d\n", params.newton_max_iter);
  printf("   lbiter:          %d\n", params.lbiter);
  printf("   lbtol:           %.5e\n", params.lbtol);
  printf("   cgmaxiter:       %d\n", params.cg_params.max_iter);
  printf("   cgtol:           %.5e\n", params.cg_params.tol);

  printf("NB: lbiter and tau usually generated automatically.\n");


  u = malloc_double(N);

  l1qc_newton(N, x, u, b, M, pix_idx, params);

  plhs[0] = mxCreateDoubleMatrix((mwSize)N, 1, mxREAL);
  x_out =  mxGetPr(plhs[0]);
  for (i=0; i<N; i++){
    x_out[i] = x[i];
  }

  free_double(u);
  free(pix_idx);

} /* ------- mexFunction ends here ----- */
