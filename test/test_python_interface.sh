#!/bin/bash

set -u
source "${srcdir}/test_runner_utils.sh"

cur_dir=$(pwd)

cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z ${verbose+x}; then
    verbose=0
fi


data_path="${TEST_DATA_DIR}/example_img_data.json"


if test x$WITH_PYTHON = xno;
then
    exit 77
fi

failures=0


python3 "${ABS_TOP_BUILDDIR}/interfaces/python/test_L1cPy.py" -v
failures+=$?

# If failures >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
