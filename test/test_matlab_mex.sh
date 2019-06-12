#!/bin/bash
set -u
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


if test -z ${verbose+x}; then
    verbose=0
fi
failures=0
# We must supply the json path as a function argument.
# The function ml_cmd() will call init_matlab_paths.m which
# will setup matlabs path for us.

data_path="${TEST_DATA_DIR}/example_img_data.json"


# Check the l1qc_dct mex interface
cmd=$(ml_cmd  "test_l1qc_dct_mex" "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
failures+=$?

# Check the breg_anistropic_TV interface
cmd=$(ml_cmd "test_breg_TV_mex" "'${data_path}'" "${verbose}")

matlab -nojvm -r "${cmd}"
failures+=$?

# If failurs >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
