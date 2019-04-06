#include "config.h"
/* It is important that mex.h is included before l1c_common.h
 so that, for this function, we get matlabs definition of fprintf
*/
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>

#include "l1c.h"
#include "l1qc_newton.h"


#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

int check_input_size(l1c_int N){
  /* Prevent segfault from bug in l1qc_newton, with splitting up DWORK.
     We can allow 512*512, but not 511*511.

     Is there a cleaner way to check this?
  */
  double divisor = (double) DALIGN / (double)sizeof(double); // e.g., 8
  double maybe_div = (double)N / divisor;

  if ((int)maybe_div == N/((int)DALIGN/(int)sizeof(double))){
    return 0;
  }

  return 1;
}

/*
 *	m e x F u n c t i o n

 The matlab protype is
 [x, LBRes] = l1qc_dct(N, M, b, pix_idx, opts),
 where
 b has length m
 pix_idx has length m
 and opts is an options struct. See the l1qc_dct.m doc string for details.
 */
void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  // l1qc(x0, b, pix_idx, params);
  /* inputs */


  double *x_out=NULL,  *b=NULL;

  l1c_L1qcOpts l1qc_opts = {.epsilon=0,
                            .mu = 0,
                            .lbtol = 0,
                            .newton_tol = 0,
                            .newton_max_iter = 0,
                            .verbose = 0,
                            .l1_tol = 0,
                            .lbiter = 0,
                            .cg_tol = 0,
                            .cg_maxiter = 0,
                            .cg_verbose = 0,
                            .warm_start_cg=0};


  l1c_LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  l1c_int i=0, M=0,Nrow=0, Ncol=1, npix=0, NM=0, idx=0;
  double *pix_idx_double=NULL, *x_ours=NULL;
  l1c_int *pix_idx;
  int status=0;
  // mwSize *dims;
  if( !(nlhs > 0)) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:nlhs",
                      "One output required.");
  }

  if( (nrhs != 5) ){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:nrhs",
                      "five inputs required.");
  }


  /* -------- Check Nrow -------------*/
  if( !mxIsDouble(prhs[0]) || !mxIsScalar(prhs[0]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "First Input vector must be type double.");
  }
  Nrow = (l1c_int) mxGetScalar(prhs[0]);


  /* -------- Check Ncol -------------*/
  if( !mxIsDouble(prhs[1]) || !mxIsScalar(prhs[1]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "First Input vector must be type double.");
  }
  Ncol = (l1c_int) mxGetScalar(prhs[1]);
  NM = Ncol * Nrow;

  /* Until I fix the DWORK/DALIGN issue, check that NM is divisible by DALIGN/sizeof(double)*/

  if (check_input_size(NM)){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:not_align-able",
                      "l1qc_dct currently requires N*M be divisible by %d.", DALIGN/sizeof(double));
  }

  /* -------- Check b -------------*/
  if( !mxIsDouble(prhs[2]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "Second Input vector must be type double.");
  }

  /* check that b is vector */
  if( (mxGetN(prhs[2]) > 1)  & (mxGetM(prhs[2]) >1) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notVector",
                      "Second Input must be a vector.");
  }else{
    M = (l1c_int) ( mxGetM(prhs[2]) *  mxGetN(prhs[2]) );
  }

  b = mxGetPr(prhs[2]);


  /* -------- Check pix_idx -------------*/
  if( !mxIsDouble(prhs[3]) ) {
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notDouble",
                      "pix_idx vector must be type double.");
  }

  /* check that pix_idx input argument is vector */
  if( (mxGetN(prhs[3]) > 1)  & (mxGetM(prhs[3]) >1) ){
    printf("num dim = %d\n", mxGetNumberOfDimensions(prhs[1]));
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notVector",
                      "Third Input must be a vector.");
  }else{
    npix = (l1c_int) ( mxGetM(prhs[3]) *  mxGetN(prhs[3]) );
  }

  pix_idx_double = mxGetPr(prhs[3]);

  if ( !(Nrow*Ncol > M) || (M != npix) || Nrow <=0 || Ncol <=0){
    printf("Nrow = %d, Ncol=%d, M=%d, npix = %d\n", Nrow, Ncol, M, npix);
    mexErrMsgIdAndTxt("l1c:l1qc_dct:incompatible_dimensions",
                      "Must have length(x0) > length(b), and length(b) = length(pix_idx).");
  }


  /* -------- Check params struct -------------*/
  if ( !mxIsStruct(prhs[4])){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:notStruct",
                      "params must be a struct.");
  }

  int nfld = mxGetNumberOfFields(prhs[4]);
  // char *flds[]={ "verbose", "tau", "mu"};
  const char *name;

  mxArray *tmp;
  for (i=0; i<nfld; i++){
    tmp = mxGetFieldByNumber(prhs[4], 0, i);
    name =mxGetFieldNameByNumber(prhs[4], i);
    if (!tmp){
      mexErrMsgIdAndTxt("l1c:l1qc_dct:norverbose",
                        "Error loading field '%s'.", name);
    }else if( !mxIsScalar(tmp) | !mxIsNumeric(tmp) ){
      mexErrMsgIdAndTxt("l1c:l1qc_dct:norverbose",
                        "Bad data in field '%s'. Fields in params struct must be numeric scalars.", name);
    }

    if ( strcmp(name, "epsilon") == 0){
      l1qc_opts.epsilon = mxGetScalar(tmp);
    }else if ( strcmp(name, "mu") == 0){
      l1qc_opts.mu = mxGetScalar(tmp);
    }else if ( strcmp(name, "tau") == 0){
      l1qc_opts.tau = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_tol") == 0){
      l1qc_opts.newton_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_max_iter") == 0){
      l1qc_opts.newton_max_iter = (int) mxGetScalar(tmp);
    }else if ( strcmp(name, "lbiter") == 0){
      l1qc_opts.lbiter = (int)mxGetScalar(tmp);
    }else if ( strcmp(name, "lbtol") == 0){
      l1qc_opts.lbtol = mxGetScalar(tmp);
    }else if ( strcmp(name, "l1_tol") == 0){
      l1qc_opts.l1_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgtol") == 0){
      l1qc_opts.cg_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgmaxiter") == 0){
      l1qc_opts.cg_maxiter = mxGetScalar(tmp);
    }else if ( strcmp(name, "warm_start_cg") == 0){
      l1qc_opts.warm_start_cg = (int) mxGetScalar(tmp);
    }else if ( strcmp(name, "verbose") == 0){
      l1qc_opts.verbose= (int) mxGetScalar(tmp);
    }else{
      mexErrMsgIdAndTxt("l1c:l1qc_dct:norverbose",
                        "Unrecognized field '%s' in params struct.", name);
    }

  }

  if (l1qc_opts.verbose > 1){
    printf("Input Parameters\n---------------\n");
    printf("   verbose:         %d\n", l1qc_opts.verbose);
    printf("   epsilon:         %.5e\n", l1qc_opts.epsilon);
    printf("   tau:             %.5e\n", l1qc_opts.tau);
    printf("   mu:              %.5e\n", l1qc_opts.mu);
    printf("   newton_tol:      %.5e\n", l1qc_opts.newton_tol);
    printf("   newton_max_iter: %d\n", l1qc_opts.newton_max_iter);
    printf("   lbiter:          %d\n", l1qc_opts.lbiter);
    printf("   lbtol:           %.5e\n", l1qc_opts.lbtol);
    printf("   l1_tol:           %.5e\n", l1qc_opts.l1_tol);
    printf("   cgmaxiter:       %d\n", l1qc_opts.cg_maxiter);
    printf("   cgtol:           %.5e\n", l1qc_opts.cg_tol);

    printf("NB: lbiter and tau usually generated automatically.\n");
  }





  /* We are going to change x, so we must allocate and make a copy, so we
     dont change data in Matlabs workspace.
  */
  x_ours = l1c_malloc_double(Nrow*Ncol);
  /*
    pix_idx will naturally be a double we supplied from matlab
    and will have 1-based indexing. Convert to integers and
    shift to 0-based indexing.
  */
  pix_idx = calloc(M, sizeof(l1c_int));
  if (!x_ours || !pix_idx){
    status = L1C_OUT_OF_MEMORY;
    goto exit;
  }
  /* We need to validate pix_idx.*/
  for (i=0; i<M; i++){
    idx = ((l1c_int) pix_idx_double[i]) - 1;
    if (idx <0 || idx > NM-1){
      status = L1C_INVALID_ARGUMENT;
      goto exit;
    }
    pix_idx[i] = idx;
  }

  int stat = l1qc_dct(Nrow, Ncol, x_ours, M, b, pix_idx, l1qc_opts, &lb_res);
  if (stat){
    lb_res.status = stat;
  }


  /* Prepare output data.
   */
  plhs[0] = mxCreateDoubleMatrix((mwSize)Nrow * (mwSize)Ncol, 1, mxREAL);

  x_out =  mxGetPr(plhs[0]);

  for (i=0; i<Nrow*Ncol; i++){
    x_out[i] = x_ours[i];
  }
  /* Only build the output struct if there is more than 1 output.*/
  if (nlhs == 2){
    const char *fnames[] = {"l1",
                            "total_newton_iter",
                            "total_cg_iter",
                            "status"};

    mxArray *l1_mex_pr, *total_newton_iter_mex_pr;
    mxArray *total_cg_iter_mex_pr, *status_mex_pr;
    plhs[1] = mxCreateStructMatrix(1, 1, NUMBER_OF_FIELDS(fnames), fnames);

    l1_mex_pr = mxCreateDoubleMatrix(1,1, mxREAL);
    total_newton_iter_mex_pr = mxCreateDoubleMatrix(1,1, mxREAL);
    total_cg_iter_mex_pr = mxCreateDoubleMatrix(1,1, mxREAL);
    status_mex_pr            = mxCreateDoubleMatrix(1,1, mxREAL);

    *mxGetPr(l1_mex_pr) = lb_res.l1;
    *mxGetPr(total_newton_iter_mex_pr) = (double)lb_res.total_newton_iter;
    *mxGetPr(total_cg_iter_mex_pr) = (double)lb_res.total_cg_iter;
    *mxGetPr(status_mex_pr) = (double)lb_res.status;

    mxSetField(plhs[1], 0, "l1", l1_mex_pr);
    mxSetField(plhs[1], 0, "total_newton_iter", total_newton_iter_mex_pr);
    mxSetField(plhs[1], 0, "total_cg_iter", total_cg_iter_mex_pr);
    mxSetField(plhs[1], 0, "status", status_mex_pr);
  }

 exit:
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

  default:
    break;
  }

} /* ------- mexFunction ends here ----- */
