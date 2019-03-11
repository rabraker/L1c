#include "config.h"
/* It is important that mex.h is included before l1c_common.h
 so that, for this function, we get matlabs definition of fprintf
*/
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <stdlib.h>
#include "l1qc_newton.h"
#include "l1c_common.h"
#include "l1c_memory.h"
#include <math.h>


/* ----- Forward Declarations */
typedef struct L1qcDctOpts{
  double epsilon;
  double mu;
  double lbtol;
  double tau;
  int lbiter;
  double newton_tol;
  int newton_max_iter;
  int verbose;
  double l1_tol;
  double cgtol;
  int cgmaxiter;
  int warm_start_cg;
}L1qcDctOpts;

// int l1qc_dct(int N, double *x_out, int M, double *b, l1c_int *pix_idx,
//              L1qcDctOpts opts, LBResult *lb_res);

int l1qc_dct(int Nrow, int Ncol, double *x_out, int M, double *b, l1c_int *pix_idx,
              L1qcDctOpts opts, LBResult *lb_res);


#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

/*
 *	m e x F u n c t i o n

 The matlab protype is
 [x, LBRes] = l1qc_dct(N, b, pix_idx, opts),
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

  L1qcDctOpts l1qc_dct_opts = {.epsilon=0,
							   .mu = 0,
							   .lbtol = 0,
							   .newton_tol = 0,
                               .newton_max_iter = 0,
                               .verbose = 0,
                               .l1_tol = 0,
                               .lbiter = 0,
                               .cgtol = 0,
                               .cgmaxiter = 0,
                               .warm_start_cg=0};


  LBResult lb_res = {.status = 0, .total_newton_iter = 0, .l1=INFINITY};

  l1c_int i=0, M=0,Nrow=0, Ncol=1, npix=0;
  double *pix_idx_double=NULL, *x_ours=NULL;
  l1c_int *pix_idx;
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
      l1qc_dct_opts.epsilon = mxGetScalar(tmp);
    }else if ( strcmp(name, "mu") == 0){
      l1qc_dct_opts.mu = mxGetScalar(tmp);
    }else if ( strcmp(name, "tau") == 0){
      l1qc_dct_opts.tau = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_tol") == 0){
      l1qc_dct_opts.newton_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "newton_max_iter") == 0){
      l1qc_dct_opts.newton_max_iter = (int) mxGetScalar(tmp);
    }else if ( strcmp(name, "lbiter") == 0){
      l1qc_dct_opts.lbiter = (int)mxGetScalar(tmp);
    }else if ( strcmp(name, "lbtol") == 0){
      l1qc_dct_opts.lbtol = mxGetScalar(tmp);
    }else if ( strcmp(name, "l1_tol") == 0){
      l1qc_dct_opts.l1_tol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgtol") == 0){
      l1qc_dct_opts.cgtol = mxGetScalar(tmp);
    }else if ( strcmp(name, "cgmaxiter") == 0){
      l1qc_dct_opts.cgmaxiter = mxGetScalar(tmp);
    }else if ( strcmp(name, "warm_start_cg") == 0){
      l1qc_dct_opts.warm_start_cg = (int) mxGetScalar(tmp);
    }else if ( strcmp(name, "verbose") == 0){
      l1qc_dct_opts.verbose= (int) mxGetScalar(tmp);
    }else{
      mexErrMsgIdAndTxt("l1c:l1qc_dct:norverbose",
                        "Unrecognized field '%s' in params struct.", name);
    }

  }

  if (l1qc_dct_opts.verbose > 1){
    printf("Input Parameters\n---------------\n");
    printf("   verbose:         %d\n", l1qc_dct_opts.verbose);
    printf("   epsilon:         %.5e\n", l1qc_dct_opts.epsilon);
    printf("   tau:             %.5e\n", l1qc_dct_opts.tau);
    printf("   mu:              %.5e\n", l1qc_dct_opts.mu);
    printf("   newton_tol:      %.5e\n", l1qc_dct_opts.newton_tol);
    printf("   newton_max_iter: %d\n", l1qc_dct_opts.newton_max_iter);
    printf("   lbiter:          %d\n", l1qc_dct_opts.lbiter);
    printf("   lbtol:           %.5e\n", l1qc_dct_opts.lbtol);
    printf("   l1_tol:           %.5e\n", l1qc_dct_opts.l1_tol);
    printf("   cgmaxiter:       %d\n", l1qc_dct_opts.cgmaxiter);
    printf("   cgtol:           %.5e\n", l1qc_dct_opts.cgtol);

    printf("NB: lbiter and tau usually generated automatically.\n");
  }





  /* We are going to change x, so we must allocate and make a copy, so we
     dont change data in Matlabs workspace.
  */
  x_ours = malloc_double(Nrow*Ncol);
  /*
    pix_idx will naturally be a double we supplied from matlab
    and will have 1-based indexing. Convert to integers and
    shift to 0-based indexing.
  */
  pix_idx = calloc(M, sizeof(l1c_int));
  if (!x_ours || !pix_idx){
    mexErrMsgIdAndTxt("l1c:l1qc_dct:outofmemory",
                      "Error Allocating memory.");
  }

  for (i=0; i<M; i++){
    pix_idx[i] = ((l1c_int) pix_idx_double[i]) - 1;
  }

  int stat = l1qc_dct(Nrow, Ncol, x_ours, M, b, pix_idx, l1qc_dct_opts, &lb_res);
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

 free(pix_idx);
 free_double(x_ours);
#if defined(_USEMKL_)
 mkl_free_buffers();
#endif
} /* ------- mexFunction ends here ----- */
