#!/bin/bash

source "${srcdir}/test_runner_utils.sh"

cur_dir=$(pwd)

function cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

host_os=$(get_os)

if test $host_os = windows || test x$WITH_MEX = xno;
then
    exit 77
fi


if test -z $verbose; then
    verbose=0
fi
failures=0
# We need to tell matlab about 3 different paths:
#  1. $lib_dir, the location of mexfiles.
#  2. $ml_script_di, this is where the driver scripts are.
#  3. $mex_src_dir, this is where, e.g., l1qc_opts.m is.
#
# Additionally, we must supply the json path as a function argument.

mex_src_dir="${ABS_TOP_SRCDIR}/interfaces/mex"
lib_dir="${ABS_TOP_BUILDDIR}/interfaces/mex"
ml_script_dir="${ABS_TOP_SRCDIR}/examples/matlab"
data_path="${TEST_DATA_DIR}/example_img_data.json"


# Check the l1qc_dct mex interface
cmd=$(ml_cmd "${mex_src_dir}" "${lib_dir}" "${ml_script_dir}" "test_l1qc_dct_mex" \
             "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
failures+=$?

# Check the breg_anistropic_TV interface
cmd=$(ml_cmd "${mex_src_dir}" "${lib_dir}" "${ml_script_dir}" "test_breg_TV_mex" \
             "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
failures+=$?

# If failurs >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
