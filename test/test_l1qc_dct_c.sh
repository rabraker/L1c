#!/bin/sh

cur_dir=$(pwd)

cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z $verbose; then
    verbose=0
fi



LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"

"${LIB_DIR}/l1qc_dct_c" "${data_path}"> /dev/null

exit $?
