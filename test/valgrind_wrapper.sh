#!/bin/bash
LT="${TOP_BUILDDIR}/libtool"
# FFTW with threads tends to leave stuff around, even after calling fftw_destroy_plan.
# However, this only shows up as possible lost in valgrind.


source "${srcdir}/test_runner_utils.sh"
host_os=$(get_os)
if test $host_os = windows;
then
    exit 77
fi


if test x$with_valgrind = xyes;
then

    VALGRIND_OPTS="--leak-check=full --errors-for-leak-kinds=definite --error-exitcode=1"
    # CK_RUN_CASE=l1qc_gradient
    $LT --mode=execute valgrind $VALGRIND_OPTS "${TOP_BUILDDIR}/test/test_l1c.test"

    data_path="${ABS_TOP_SRCDIR}/test/test_data/example_img_data127.json"
    $LT --mode=execute valgrind $VALGRIND_OPTS "${TOP_BUILDDIR}/examples/c/l1qc_dct_c" "${data_path}"
else
    exit 77
fi
