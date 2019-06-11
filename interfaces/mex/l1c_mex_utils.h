#ifndef __L1C_MEX_UTILS__
#define __L1C_MEX_UTILS__

#include "mex.h"
#include "matrix.h"
#include "l1c.h"

#define NUMBER_OF_FIELDS(ST) (sizeof(ST)/sizeof(*ST))

int
check_input_size(l1c_int N);

void
_mex_assert_input_double(const mxArray ** mxa, size_t idx);

void
_mex_assert_num_outputs(int N, int N_required);

void
_mex_assert_num_inputs(int N, int N_required);

double
_mex_get_double_scaler_or_fail(const mxArray **mxa, size_t idx);

double
_mex_get_int32_scaler_or_fail(const mxArray **mxa, size_t idx);

void
_mex_assert_1darray(const mxArray **mxa, size_t idx);

void
_mex_assert_2Darray_with_size(const mxArray **mxa, size_t idx, size_t n, size_t m);

void
_mex_get_double_array_or_fail(const mxArray **mxa, size_t idx, double *x[], size_t *N);

void
_mex_get_int_array_or_fail(const mxArray **mxa, size_t idx, int *x[], size_t *N);

#endif
