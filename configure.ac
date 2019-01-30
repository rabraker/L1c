#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# autoconf && CPPFLAGS=-I${MKLROOT}/include LDFLAGS=-L${MKLROOT}/lib/intel64 ./configure --with-mkl --with-fftw3
AC_PREREQ([2.69])
AC_INIT([L1c], [0.1], [abraker@rabraker.com])

AC_CONFIG_SRCDIR([src/vcl/vcl_math.cpp])
AC_CONFIG_SRCDIR([src/l1c_math.c])

dnl We do not want -g -O2 to be the default.
dnl But dont over-ride the user. These must
dnl come before AC_PROG_CXX etc.

if test -z $CFLAGS;
then
    CFLAGS=''
fi
if test -z $CXXLAGS;
then
    CXXFLAGS=''
fi

AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

dnl AC_PROG_MAKE_SET
dnl AC_PROG_CPP
AM_PROG_AR
dnl AM_PROG_LIBTOOL
dnl AC_PROG_INSTALL

AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
dnl This is necessary, otherwise, when we use AC_CHECK_HEADER inside the conditional,
dnl AC_PROG_CPP will be expanded in the else, and not visible outside, which leads to
dnl warnings like WARNING: check.h: accepted by the compiler, rejected by the preprocessor!
AC_PROG_CPP



AC_ARG_WITH(debug, [AS_HELP_STRING([--with-debug],[Compile with debugging support.])],
                   [AC_DEFINE(WITH_DEBUG,1,0)
                   AC_SUBST(DBG, "-O0 -ggdb")],
                   [])


AC_ARG_WITH(opt, [AS_HELP_STRING([--without-opt],[Compile optimizations off.])],
[OPT=no], [OPT=yes])

# ------------------------------ Checks for MKL -------------------------------------------
AC_ARG_WITH(mkl, [AS_HELP_STRING([--with-mkl], [Compile using Intel Math Kernel Library])], [
                 use_mkl=yes], [use_mkl=no])


echo $use_mkl
if test "x$use_mkl" = "xyes";
then
    AX_CHECK_LINK_FLAG([-lmkl_rt], [AX_APPEND_FLAG([-lmkl_rt], [L1C_LDLIBS])],
                                   [echo "Cannot find MKL"],[], [])

    AX_APPEND_ALL_COMPILE_FLAGS([-fopenmp -m64 -DMKL_ILP64 -D_USEMKL_], [TMP],
                                   [AC_MSG_ERROR("Compiler does not support required flag")],
                                   [],[])
    AC_SUBST(L1CCPPFLAGS, $TMP)
    AC_CHECK_HEADER([mkl.h], [], [AC_MSG_ERROR("Cannot find mkl.h")])
    dnl # FIXME: Replace `main' with a function in `-lpthread':
    AC_CHECK_LIB([pthread], [pthread_mutex_trylock], [],
    AC_MSG_ERROR(Cannot find libpthread))
    echo stuff
else
 dnl look for another BLAS library
     AC_CHECK_LIB([openblas], [cblas_daxpby],[],[
                 AC_MSG_ERROR("BLAS library with cblas_daxpby is required for this program")])

fi



echo "Compile with MKL? $use_mkl"
echo $CFLAGS $CPPFLAGS




# ------------------------------ Checks for FFTW ---------------------------------------- #
AC_ARG_WITH(fftw3, [AS_HELP_STRING([--with-fftw3], [Compile example with FFTW])], [
                 use_fftw=yes], [use_fftw=no])

if test x$use_fftw = xyes;
then
    AX_CHECK_LINK_FLAG([-lfftw3], [AX_APPEND_FLAG([-lmkl_rt], [L1C_LDLIBS])],
                                  [echo "Cannot find MKL"],[], [])

   AC_CHECK_HEADER([fftw3.h], [], [AC_MSG_ERROR("Cannot find fftw3.h")])
   AC_CHECK_LIB([fftw3], [fftw_plan_r2r_1d],[],[
                echo "fftw3 library is required for this program"
                exit -1])
fi

echo "Compile with fftw3? $use_fftw"


AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wextra -Wunused -Werror  \
                         -Wcast-qual -Wcast-align -Winit-self      \
                         -Wfloat-equal -Wundef -Wshadow -Wswitch-enum -Wvla], [WARNFLAGS], [] [])
AC_SUBST(WARNFLAGS, $WARNFLAGS)
echo  $CFLAGS $CPPFLAGS


dnl echo "WARN_FLGS = $WARN_FLAGS"
opt_enable=yes
if test x$OPT = xyes ;
then
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [OPT], [], [])
fi
echo "OPT = $OPT"

# Checks for programs.


# Checks for libraries.
AC_CHECK_LIB([check], [suite_add_tcase],[],[
             echo "check library is required for this program"
             exit -1])
AC_CHECK_LIB([cjson], [cJSON_GetArraySize],[],[
    echo "cJSON library is required for this program"
    exit -1])


AC_CHECK_LIB([gomp], [omp_set_num_threads],[],[
echo "omp library is required for this program"
exit -1])

AC_CHECK_LIB([m], [sqrt],[],[
   echo "math library is required for this program"
   exit -1])



# Check for headers

AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])
AC_CHECK_HEADER([check.h], [], [AC_MSG_ERROR("Cannot find check.h")], [])
AC_CHECK_HEADER([check.h], [], [AC_MSG_ERROR("Cannot find check.h")])
AC_CHECK_HEADER([cjson/cJSON.h], [], [AC_MSG_ERROR("Cannot find CJSON.h")])
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR("Cannot find math.h")])


dnl # FIXME: Replace `main' with a function in `-ldl':
dnl AC_CHECK_LIB([dl], [main])




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


AC_CONFIG_FILES([Makefile  src/Makefile src/vcl/Makefile interfaces/Makefile])
AC_OUTPUT
