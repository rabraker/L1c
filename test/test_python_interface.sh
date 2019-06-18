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

#N.B. Somehow, by almost pure magic, when we read these environmental
# variables in python in mingw, the path will be converted to a windows
# style path, e.g., C:/msys64/...
export L1C_INTERFACE_DIR="${ABS_TOP_BUILDDIR}/interfaces/python"
export L1C_SRC_LIB_DIR="${ABS_TOP_BUILDDIR}/src/.libs"

python3 "${ABS_TOP_BUILDDIR}/interfaces/python/examples/dct_example.py" "${data_path}"

failures+=$?


python3 "${ABS_TOP_BUILDDIR}/interfaces/python/test_L1cPy.py" -v
failures+=$?

# If failures >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
