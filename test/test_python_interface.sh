#!/bin/bash

source "${srcdir}/test_runner_utils.sh"

cur_dir=$(pwd)

cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z $verbose; then
    verbose=0
fi


data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"


# host_os=$(get_os)
# if test $host_os = windows;
# then
#     exit 77
# fi



#N.B. Somehow, by almost pure magic, when we read these environmental
# variables in python in mingw, the path will be converted to a windows
# style path, e.g., C:/msys64/...
export L1C_INTERFACE_DIR="${ABS_TOP_BUILDDIR}/interfaces"
export L1C_SRC_LIB_DIR="${ABS_TOP_BUILDDIR}/src/.libs"

python3 "${ABS_TOP_SRCDIR}/interfaces/dct_example.py" "${data_path}"

exit $?
