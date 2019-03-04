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
    $LT --mode=execute valgrind $VALGRIND_OPTS ./test_l1c
else
    exit 77
fi
