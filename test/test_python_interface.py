#!/bin/sh

source "${srcdir}/test_runner_utils.sh"

cur_dir=$(pwd)

function cleanup(){
    cd $cur_dir
}
trap cleanup EXIT

if test -z $verbose; then
    verbose=2
fi



LIB_DIR="${ABS_TOP_BUILDDIR}/interfaces"
data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data.json"

# host_os=$(get_os)
# if test $host_os = windows;
# then
#     exit 77
# fi


# # check the python dlopen interface
LIB_DIR="${LIB_DIR}" python3 "${ABS_TOP_SRCDIR}/interfaces/dct_example.py" "${data_path}" 0

exit $?
