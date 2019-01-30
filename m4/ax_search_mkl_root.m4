# ============================================================================
#  https://www.gnu.org/software/autoconf-archive/ax_append_compile_flags.html
# ============================================================================
#
# SYNOPSIS
#
#   AX_APPEND_COMPILE_FLAGS([FLAG1 FLAG2 ...], [FLAGS-VARIABLE], [EXTRA-FLAGS], [INPUT])
#
# DESCRIPTION
#
#   Search for Intel Math Kernel Library in standard locations.
#
# LICENSE
#
#   Copyright (c) 2011 Maarten Bosmans <mkbosmans@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 7

AC_DEFUN([AX_SEARCH_MKL_ROOT],[dnl

# We start simple: The default install location seems to be /opt/intel
echo $1
if test x$1 = xyes;
then
echo "Trying to find MKL"
     MKL_ROOT=/opt/intel/mkl/
elif test ! -z "$1";
then
     MKL_ROOT="$1"
     echo "maybe: $maybe_mkl_root"
fi


echo $MKL_ROOT
if test ! -e ${MKL_ROOT}/include/mkl.h;
then
    AC_MSG_ERROR("Could not locate MKL")
fi

echo "found: $MKL_ROOT\nChecking libmkl_rt"

# Save LDFLAGS and restore
TMP_LDFLAGS=$LDFLAGS
LDFLAGS=-L$MKL_ROOT/lib/intel64

AC_CHECK_LIB(mkl_rt, d_forward_trig_transform, [],
            [AC_MSG_ERROR("Could not locate MKL")])

LDFLAGS=$TMP_LDFLAGS


AM_CPPFLAGS=-I${MKL_ROOT}/include
AC_SUBST(MKL_LDFLAGS, $-L{MKL_ROOT}/lib/intel64)
AC_SUBST(MKL_LDLIBS, -lmkl_rt)




])dnl AX_SEARCH_MKL_ROOT
