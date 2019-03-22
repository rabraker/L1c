#!/bin/bash

cur_dir=$(pwd)

cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z $verbose; then
    verbose=0
fi


failures=0

LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"

data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"
"${LIB_DIR}/l1qc_dct_c" "${data_path}"

failures+=$?

# Regression check that we dont have to have N divisible by
# DALIGN/sizeof(double = 64/8
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data127.json"
"${LIB_DIR}/l1qc_dct_c" "${data_path}"

failures+=$?

# If failures >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
