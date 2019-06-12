#include <stdint.h>
#include "mex.h"
#include "matrix.h"
#include "l1c_mex_utils.h"



int
check_input_size(l1c_int N){
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

void
_mex_assert_input_double(const mxArray ** mxa, size_t idx){
  if( !mxIsDouble(mxa[idx]) ) {
    mexErrMsgIdAndTxt("l1c:notDouble",
                      "Argument %d must be a double.", idx+1);
  }
}


void
_mex_assert_input_double_scalar(const mxArray ** mxa, size_t idx){
  if( !mxIsDouble(mxa[idx]) || !mxIsScalar(mxa[idx]) ) {
    mexErrMsgIdAndTxt("l1c:notDoubleScaler",
                      "Argument %d must be a scaler double.", idx+1);
  }
}

void
_mex_assert_num_outputs(int N, int N_required){
  if(N < N_required ) {
    mexErrMsgIdAndTxt("l1c:nlhs", "%d output(s) required.", N_required);
  }
}


void
_mex_assert_num_inputs(int N, int N_required){
    if( N < N_required ) {
      mexErrMsgIdAndTxt("l1c:nlhs", "%d input(s)  required.", N_required);
    }
}


double
_mex_get_double_scaler_or_fail(const mxArray **mxa, size_t idx){
  _mex_assert_input_double_scalar(mxa, idx);

  return (double) mxGetScalar(mxa[idx]);
}

l1c_int
_mex_get_int32_scaler_or_fail(const mxArray **mxa, size_t idx){
  if( !mxIsInt32(mxa[idx]) || !mxIsScalar(mxa[idx]) ) {
    mexErrMsgIdAndTxt("l1c:notInt32Scaler",
                      "Argument %d must be a scaler double.", idx+1);
  }
  return (l1c_int) mxGetScalar(mxa[idx]);
}

void
_mex_assert_1darray(const mxArray **mxa, size_t idx){
    size_t n = mxGetN(mxa[idx]);
    size_t m = mxGetM(mxa[idx]);
    size_t d = mxGetNumberOfDimensions(mxa[idx]);

    if ( (n==0) || (m==0) ){
       mexErrMsgIdAndTxt("l1c:not1DArray", "Input %d must not be empty.", idx+1);
    }
    if (mxIsScalar(mxa[idx])){
      mexErrMsgIdAndTxt("l1c:not1DArray", "Input %d must not be a scalar.", idx+1);
    }
    /* In matlab, [1,2,3], or [1;2;3], will both have d==2. */
    if(d > 2 || ( (n!=1) && (m!=1))){
      mexErrMsgIdAndTxt("l1c:not1DArray", "Input %d must have only 1 dimension.", idx+1);
    }
}


void
_mex_assert_2Darray_with_size(const mxArray **mxa, size_t idx, l1c_int n, l1c_int m){
  /* Matlab stores as column major. But the call M=number of rows and N=number of columns
   for some bizzar reason. */
  size_t n_act = mxGetM(mxa[idx]);
  size_t m_act = mxGetN(mxa[idx]);
  size_t n_sizet = (size_t) n;
  size_t m_sizet = (size_t) m;
  size_t d = mxGetNumberOfDimensions(mxa[idx]);

  if ( (n_act==0) || (m_act==0) ){
    mexErrMsgIdAndTxt("l1c:not2DArray", "Input %d must not be empty.", idx+1);
  }
  if (mxIsScalar(mxa[idx])){
    mexErrMsgIdAndTxt("l1c:not2DArray", "Input %d must not be a scalar.", idx+1);
  }

  if(d > 2 ){
    mexErrMsgIdAndTxt("l1c:not2DArray", "Input %d must have only 1 dimension.", idx+1);
  }
  if (n_act < n_sizet || m_act < m_sizet){
    mexErrMsgIdAndTxt("l1c:ArraySize", "Input %d has size (%d, %d); expected at least (%d, %d)\n",
                      idx+1, n_act, m_act, n, m);
  }
}


void
_mex_get_double_array_or_fail(const mxArray **mxa, size_t idx, double *x[], l1c_int *N){
  size_t n = mxGetN(mxa[idx]);
  size_t m = mxGetM(mxa[idx]);
  size_t nm=0;

  if (!mxIsDouble(mxa[idx])){
    mexErrMsgIdAndTxt("l1c:notDouble", "Input %d must be a double.", idx);
  }
  _mex_assert_1darray(mxa, idx);

  nm = n * m;
  if (nm > INT_MAX){
    mexErrMsgIdAndTxt("l1c:ArrayTooLarge",
                      "L1c can only index up to %d\n", INT_MAX);
  }

  *N = (l1c_int) nm;
  *x = mxGetPr(mxa[idx]);
}


void
_mex_get_int_array_or_fail(const mxArray **mxa, size_t idx, int *x[], l1c_int *N){
  size_t n = mxGetN(mxa[idx]);
  size_t m = mxGetM(mxa[idx]);
  size_t nm=0;

  if (!mxIsInt32(mxa[idx])){
    mexErrMsgIdAndTxt("l1c:notInt32", "Input %d must be an int32.", idx);
  }
  _mex_assert_1darray(mxa, idx);

  nm = n * m;
  if (nm > INT_MAX){
    mexErrMsgIdAndTxt("l1c:ArrayTooLarge",
                      "L1c can only index up to %d\n", INT_MAX);
  }

  *N = (l1c_int) nm;
  *x = mxGetInt32s(mxa[idx]);

}

