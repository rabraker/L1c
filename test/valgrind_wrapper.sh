#!/bin/bash
LT="${ABS_TOP_BUILDDIR}/libtool"
# FFTW with threads tends to leave stuff around, even after calling fftw_destroy_plan.
# However, this only shows up as "possibly lost" in valgrind.


source "${srcdir}/test_runner_utils.sh"
host_os=$(get_os)
if test $host_os = windows;
then
    exit 77
fi


if test x$with_valgrind = xyes;
then

    VALGRIND_OPTS="--leak-check=full --errors-for-leak-kinds=definite --error-exitcode=1"
    $LT --mode=execute valgrind $VALGRIND_OPTS "${ABS_TOP_BUILDDIR}/test/test_l1c.test"

    data_path="${TEST_DATA_DIR}/example_img_data127.json"
    $LT --mode=execute valgrind $VALGRIND_OPTS "${ABS_TOP_BUILDDIR}/examples/c/l1qc_dct_c" "${data_path}"
else
    exit 77
fi
