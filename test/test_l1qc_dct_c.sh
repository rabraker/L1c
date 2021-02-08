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

if [[ x$ABS_TOP_BUILDDIR = x ]]; then
    echo "Environment variable ABS_TOP_BUILDDIR not set"
    exit 1
fi
if [[ x$TEST_DATA_DIR = x ]]; then
    echo "Environment variable TEST_DATA_DIR not set"
    exit 1
fi


BIN_DIR="${ABS_TOP_BUILDDIR}/examples/c"

data_path="${TEST_DATA_DIR}/example_img_data.json"
"${BIN_DIR}/l1qc_dct_c" "${data_path}" $verbose

failures+=$?

# Regression check that we dont have to have N divisible by
# DALIGN/sizeof(double = 64/8
data_path="${TEST_DATA_DIR}/example_img_data127.json"
"${BIN_DIR}/l1qc_dct_c" "${data_path}" $verbose

failures+=$?

# If failures >255, automake believes its zero, evidently.
if test $failures -gt 0; then
    exit 1
else
    exit $failures
fi
