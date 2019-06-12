#!/bin/bash

source "${ABS_TOP_SRCDIR}/test/test_runner_utils.sh"

mex_test_src_path="${ABS_TOP_SRCDIR}/interfaces/mex/test"
mex_test_build_path="${ABS_TOP_BUILDDIR}/interfaces/mex/test"

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
mat_cmd=$(cat <<EOF
try; \
    addpath("${mex_test_build_path}"); \
    addpath("${mex_test_src_path}"); \
    test_mex_utils_runner; \
    exit(status); \
catch ME; \
    fprintf('%s\n', ME.message'); \
    exit(1); \
end
EOF
)

echo $mat_cmd


matlab -nojvm -r "$mat_cmd"
