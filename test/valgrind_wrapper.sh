#!/bin/bash

SCRIPT_DIR=$(dirname $(readlink -f "$0"))
TEST_DATA_DIR="$SCRIPT_DIR/test_data"

echo "script dir: $SCRIPT_DIR"

source "${SCRIPT_DIR}/test_runner_utils.sh"
host_os=$(get_os)
if test $host_os = windows;
then
    echo "Windows does not support valgrind."
    exit 77
fi

# FFTW with threads tends to leave stuff around, even after calling fftw_destroy_plan.
# However, this only shows up as "possibly lost" in valgrind.
VALGRIND_OPTS="--leak-check=full --errors-for-leak-kinds=definite --error-exitcode=1"
data_path="${TEST_DATA_DIR}/example_img_data127.json"
echo "data_path= ${data_path}"

if [ x$ABS_TOP_BUILD_DIR != x ]; then
   # We are running with autotools.
   LT="${ABS_TOP_BUILDDIR}/libtool --mode=execute"
else
    LT=""
fi


if test x$with_valgrind = xyes; then
    export TEST_DATA_DIR=$TEST_DATA_DIR
    $LT valgrind $VALGRIND_OPTS "${SCRIPT_DIR}/test_l1c.test"

    $LT valgrind $VALGRIND_OPTS "${SCRIPT_DIR}/../examples/c/l1qc_dct_c" "${data_path}"
else
    echo "with_valgrind not set, skipping this test."
    exit 77
fi
