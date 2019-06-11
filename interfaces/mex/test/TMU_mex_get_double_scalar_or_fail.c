#include "config.h"
#include "mex.h"
#include "matrix.h"

#include "l1c_mex_utils.h"


void  mexFunction( int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[] )
{
  (void) nlhs;
  (void) plhs;
  (void) prhs;
  _mex_assert_num_inputs(nrhs, 1);

  double res = _mex_get_double_scaler_or_fail(prhs, 1);
  (void) res;
}
