#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include "l1qc_newton.h"
#include "l1c_common.h"
#include <math.h>

#include "dct.h" /* Or whatever you want to use for
                    the Ax Aty transformations */

#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

/*
 *	m e x F u n c t i o n
 */
void mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{

#ifdef _USEMKL_
  /* Ensure intel doesnt fuck us.*/
  mkl_set_interface_layer(MKL_INTERFACE_ILP64);
  mkl_set_threading_layer(MKL_THREADING_GNU);
#endif
  // l1qc(x0, b, pix_idx, params);
  /* inputs */
  AxFuns Ax_funs = {.Ax=dct_EMx_new,
                    .Aty=dct_MtEty,
                    .AtAx=dct_MtEt_EMx_new};

  double *x_theirs=NULL, *x_ours=NULL, *u=NULL, *b=NULL;
  NewtParams params = {.epsilon=0, .tau=0, .mu=0,
                        .newton_tol=0, .newton_max_iter = 0, .lbiter=0,
                        .lbtol=0, .verbose = 0, .cg_params.tol=0};
  LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  int i=0,N=0, M=0, npix=0;
  double *pix_idx_double=NULL, *x_out=NULL;
  l1c_int *pix_idx;
  // mwSize *dims;

  if(nrhs != 4) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:nrhs",
                      "four inputs required.");
  }

  if( !(nlhs > 0)) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:nlhs",
                      "One output required.");
  }

  /* -------- Check x0 -------------*/
  if( !mxIsDouble(prhs[0]) ) {
      mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "First Input vector must be type double.");
  }

  /* check that number of rows in second input argument is 1 */
  if( (mxGetN(prhs[0]) > 1)  & (mxGetM(prhs[0]) >1) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notVector",
                      "First Input must be a vector.");
  }else{
    N = mxGetM(prhs[0]) *  mxGetN(prhs[0]);
  }

  /* -------- Check b -------------*/
  if( !mxIsDouble(prhs[1]) ) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "Second Input vector must be type double.");
  }

  /* check that second input argument is vector */
  if( (mxGetN(prhs[1]) > 1)  & (mxGetM(prhs[1]) >1) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notVector",
                      "Second Input must be a vector.");
  }else{
    M = mxGetM(prhs[1]) *  mxGetN(prhs[1]);
  }


  /* -------- Check pix_idx -------------*/
  if( !mxIsDouble(prhs[2]) ) {
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notDouble",
                      "pix_idx vector must be type double.");
  }

  /* check that pix_idx input argument is vector */
  if( (mxGetN(prhs[2]) > 1)  & (mxGetM(prhs[2]) >1) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notVector",
                      "Third Input must be a vector.");
  }else{
    npix = mxGetM(prhs[2]) *  mxGetN(prhs[2]);
  }


 if ( !(N > M) | (M != npix) ){
   printf("N = %d, M=%d, npix = %d\n", N, M, npix);
   mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:incompatible_dimensions",
                     "Must have length(x0) > length(b), and length(b) = length(pix_idx).");
 }

 if ( !mxIsStruct(prhs[3])){
   mexErrMsgIdAndTxt("l1qc:l1qc_log_barrier:notStruct",
                     "params must be a struct.");
 }


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

  if (params.verbose > 1){
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
  }

  x_theirs = mxGetPr(prhs[0]);
  b = mxGetPr(prhs[1]);
  pix_idx_double = mxGetPr(prhs[2]);


  x_ours = malloc_double(N);
  u = malloc_double(N);

  for (i=0; i<N; i++){
    x_ours[i] = x_theirs[i];
  }

  pix_idx = calloc(M, sizeof(l1c_int));
  for (i=0; i<M; i++){
    pix_idx[i] = (l1c_int) pix_idx_double[i];
  }

  /* ---------------------------------------  */
  dct_setup(N, M, pix_idx);
  lb_res = l1qc_newton(N, x_ours, M, b,  params, Ax_funs);
  dct_destroy();

  plhs[0] = mxCreateDoubleMatrix((mwSize)N, 1, mxREAL);



  x_out =  mxGetPr(plhs[0]);

  for (i=0; i<N; i++){
    x_out[i] = x_ours[i];
  }

  if (nlhs == 2){
    const char *fnames[] = {"l1", "total_newton_iter", "status"};
    mxArray *l1_mex_pr, *total_newton_iter_mex_pr, *status_mex_pr;
    plhs[1] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS(fnames), fnames);

    l1_mex_pr = mxCreateDoubleMatrix(1,1, mxREAL);
    total_newton_iter_mex_pr = mxCreateDoubleMatrix(1,1, mxREAL);
    status_mex_pr            = mxCreateDoubleMatrix(1,1, mxREAL);

    *mxGetPr(l1_mex_pr) = lb_res.l1;
    *mxGetPr(total_newton_iter_mex_pr) = (double)lb_res.total_newton_iter;
    *mxGetPr(status_mex_pr) = (double)lb_res.status;

    mxSetField(plhs[1], 0, "l1", l1_mex_pr);
    mxSetField(plhs[1], 0, "total_newton_iter", total_newton_iter_mex_pr);
    mxSetField(plhs[1], 0, "status", status_mex_pr);
  }

  free_double(u);
  free(pix_idx);

#ifdef _USEMKL_
  mkl_free_buffers();
#endif
} /* ------- mexFunction ends here ----- */
