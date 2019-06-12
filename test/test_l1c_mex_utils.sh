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

mat_cmd=$(cat <<EOF
try; \
    init_matlab_paths;    \
    test_mex_utils; \
    exit(status); \
catch ME; \
    fprintf('%s\n', ME.message'); \
    exit(1); \
end
EOF
)

matlab -nojvm -r "$mat_cmd"
