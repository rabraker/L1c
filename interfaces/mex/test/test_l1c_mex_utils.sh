#!/bin/bash

source "${ABS_TOP_SRCDIR}/test/test_runner_utils.sh"

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

matlab -nojvm -r "try; test_mex_utils_runner; exit(0); catch; exit(1); end;"
