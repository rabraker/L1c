#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([L1c], [0.1], [abraker@rabraker.com])

AC_CONFIG_SRCDIR([src/vcl/vcl_math.cpp])
AC_CONFIG_SRCDIR([src/l1c_math.c])

dnl We do not want -g -O2 to be the default.
dnl These must come before AC_PROG_CXX etc.

CFLAGS=''
CPPFLAGS=''
CXXFLAGS=''

AC_ARG_WITH(debug, [AS_HELP_STRING([--with-debug],[add the debugging module])],
                   [AC_DEFINE(WITH_DEBUG,1,0)
                   AC_SUBST(WITH_DEBUG,1)
                   AC_SUBST(OPT,"")
                   AC_SUBST(DBG, "-O0 -ggdb")],
                   [AC_SUBST(DBG,"")])




AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

dnl AC_PROG_MAKE_SET
dnl AC_PROG_CPP

AM_PROG_AR
dnl AM_PROG_LIBTOOL
dnl AC_PROG_INSTALL


AC_SUBST(MKLROOT, /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl)

AC_SUBST(OPT," -O3 -march=native -mtune=native")
AC_SUBST(MKL_INCLUDE, "-I${MKLROOT}/include")
AC_SUBST(MKL_CFLAGS, "-DMKL_ILP64  -fopenmp -D_USEMKL_  -m64")

dnl INCLUDES ="-Iinclude -I${MKLROOT}/include"

# Checks for libraries.
dnl AC_CHECK_LIB([check], [suite_add_tcase],[],[
dnl             echo "check library is required for this program"
dnl             exit -1])

dnl AC_CHECK_LIB([cjson], [cJSON_GetArraySize],[],[
dnl     echo "cJSON library is required for this program"
dnl     exit -1])

dnl # FIXME: Replace `main' with a function in `-ldl':
dnl AC_CHECK_LIB([dl], [main])

dnl AC_CHECK_LIB([fftw3], [fftw_plan_r2r_1d],[],[
dnl    echo "fftw3 library is required for this program"
dnl    exit -1])

dnl AC_CHECK_LIB([gomp], [omp_set_num_threads],[],[
dnl    echo "omp library is required for this program"
dnl    exit -1])

dnl AC_CHECK_LIB([m], [sqrt],[],[
dnl    echo "math library is required for this program"
dnl    exit -1])

dnl AC_CHECK_LIB([mat], [matOpen],[],[
dnl    echo "libmat (Mathworks) library is required for this program"
dnl    exit -1])

dnl AC_CHECK_LIB([mex], [ mexErrMsgTxt],[],[
dnl    echo "Matlab libmex library not found."
dnl    exit -1])

# FIXME: Replace `main' with a function in `-lmwbuiltinsutil':
#AC_CHECK_LIB([mwbuiltinsutil], [main])
# FIXME: Replace `main' with a function in `-lmwservices':
#AC_CHECK_LIB([mwservices], [main])
# FIXME: Replace `main' with a function in `-lmx':
#AC_CHECK_LIB([mx], [main])

dnl AC_CHECK_LIB([openblas], [cblas_daxpby],[],[
dnl echo "A BLAS library with cblas_daxpby is required for this program"
dnl exit -1])

dnl # FIXME: Replace `main' with a function in `-lpthread':
dnl AC_CHECK_LIB([pthread], [main])

# Checks for header files.
dnl AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INT8_T
AC_C_RESTRICT
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_CHECK_FUNCS([strerror])


EXEEXT=.mexa64;

AC_CONFIG_FILES([Makefile  src/Makefile src/vcl/Makefile interfaces/Makefile])
AC_OUTPUT
