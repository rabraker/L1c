#!/bin/bash
echo $TOP_BUILDDIR
echo $(pwd)
echo "$srcdir"
LT="${TOP_BUILDDIR}/libtool"
echo "------------------ ${LT}"
# FFTW with threads tends to leave stuff around, even after calling fftw_destroy_plan.
# However, this only shows up as possible lost in valgrind.

if test x$with_valgrind = xyes;
then

    VALGRIND_OPTS="--leak-check=full --errors-for-leak-kinds=definite --error-exitcode=1"
    # CK_RUN_CASE=l1qc_gradient
    $LT --mode=execute valgrind $VALGRIND_OPTS ./test_l1c
fi
