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

LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"
TEST_SRCDIR="${ABS_TOP_SRCDIR}/test"
ml_script_dir="${TEST_SRCDIR}"
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"


# Check the l1qc_dct mex interface
cmd=$(ml_cmd "${LIB_DIR}" "${ml_script_dir}" "test_l1qc_dct_mex" \
             "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
failures+=$?

# Check the breg_anistropic_TV interface
cmd=$(ml_cmd "${LIB_DIR}" "${ml_script_dir}" "test_breg_TV_mex" \
             "'${data_path}'" "${verbose}")
matlab -nojvm -r "${cmd}"
failures+=$?

exit $failures
