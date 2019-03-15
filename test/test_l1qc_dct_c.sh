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
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json" >/dev/null

"${LIB_DIR}/l1qc_dct_c" "${data_path}"

failures+=$?

LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data127.json"

"${LIB_DIR}/l1qc_dct_c" "${data_path}"
# Right now, we expect this to fail.

failures=$?
# echo "stat = $stat"
# if test $stat -eq 0; then
#     failures+=1
# fi
# echo "---- $failures"
exit $failures
