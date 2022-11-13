#include "config.h"
#include "mex.h"
#ifdef L1C_MEX_MATLAB
#include "matrix.h"
#endif

#include "l1c_mex_utils.h"

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  (void)nlhs;
  (void)plhs;
  (void)prhs;
  _mex_assert_num_inputs(nrhs, 2);
}
