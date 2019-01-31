#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# autoconf && CPPFLAGS=-I${MKLROOT}/include LDFLAGS=-L${MKLROOT}/lib/intel64 ./configure --with-mkl --with-fftw3
#
# CPPFLAGS=-I/usr/local/OpenBlas/include LDFLAGS=-L/usr/local/OpenBlas/lib ./configure --verbose --with-fftw3
# ---- Override help
# https://stackoverflow.com/questions/29928799/can-you-override-the-default-configure-help-message-generated-as-part-of-running
dnl
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
dnl Need this for AC_EGREP_HEADER to work
AC_PROG_EGREP
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
[with_opt=yes], [with_opt=yes])


AC_ARG_WITH(mkl, [AS_HELP_STRING([--with-mkl], [Compile using Intel Math Kernel Library])], [
                 use_mkl=yes], [use_mkl=no])

AC_ARG_WITH(mex, [AS_HELP_STRING([--with-mex],
                 [Compile matlab interface mex functions. This requires that
                 CPPLAGS contains -I${MATLABROOT}/extern/include and that
                 LDFLAGS contains -L${MATLABROOT}/bin/(glnxa64|win64)/])], [
                 use_mex=yes], [use_mex=no])


# ----------------------------------------------------------------------------------------- #
# ------------------------------ Checks for MKL ------------------------------------------- #
if test "x$use_mkl" = "xyes";
then
    AX_CHECK_LINK_FLAG([-lmkl_rt], [AX_APPEND_FLAG([-lmkl_rt], [MATH_LDLIBS])],
                                   [AC_MSG_ERROR(Cannot find MKL)],[], [])

    AX_APPEND_ALL_COMPILE_FLAGS([-fopenmp -m64], [TMP],
                                   [AC_MSG_ERROR("Compiler does not support required flag")],
                                   [],[])
    AC_SUBST(L1CCPPFLAGS, $TMP)
    AC_CHECK_HEADER([mkl.h], [], [AC_MSG_ERROR("Cannot find mkl.h")])
    dnl # FIXME: Replace `main' with a function in `-lpthread':
    AC_CHECK_LIB([pthread], [pthread_mutex_trylock], [],
    AC_MSG_ERROR(Cannot find libpthread))

    AC_DEFINE([_USEMKL_], [1], 0)
    AC_DEFINE([MKL_ILP64], [1], 0)

else
 dnl ----------- look for another BLAS library ---------------------
    AC_CHECK_LIB([openblas], [cblas_daxpby],[],[
                 AC_MSG_ERROR("BLAS library with cblas_daxpby is required for this program")])

    AC_CHECK_HEADER([cblas.h], [], [AC_MSG_ERROR("Cannot find cblas.h")])
    AC_MSG_CHECKING(checking for cblas_daxpby in cblas.h)
    AC_EGREP_HEADER([cblas_daxpby], [cblas.h], [AC_MSG_RESULT(yes)],
                  [AC_MSG_RESULT(no)
                  AC_MSG_ERROR(Could not find cblas_daxpy in cblas.h)])

    AX_APPEND_FLAG([-lopenblas], [MATH_LDLIBS])

fi



echo $CFLAGS $CPPFLAGS




# ------------------------------ Checks for FFTW ---------------------------------------- #
AC_ARG_WITH(fftw3, [AS_HELP_STRING([--with-fftw3], [Compile example with FFTW])], [
                 use_fftw3=yes], [use_fftw3=no])

if test x$use_fftw3 = xyes;
then
    AX_CHECK_LINK_FLAG([-lfftw3], [AX_APPEND_FLAG([-lmkl_rt], [MATH_LDLIBS])],
                       [echo "Cannot find MKL"],[], [])

   AC_CHECK_HEADER([fftw3.h], [], [AC_MSG_ERROR("Cannot find fftw3.h")])
   AC_CHECK_LIB([fftw3], [fftw_plan_r2r_1d],[],[
                AC_MSG_ERROR(fftw3 library is required for this program)])
   AX_APPEND_FLAG([-lfftw3], [MATH_LDLIBS])
   AC_DEFINE([_USEFFTW3_], [1], 0)
fi

echo "Compile with fftw3? $use_fftw3"
AC_SUBST(MATHLIBS, $MATH_LDLIBS)


# ----------------------------------------------------------------------------------------- #
# ------------------------------ Checks for matlab MEX------------------------------------- #

if test x$use_mex = xyes;
then
dnl AC_CHECK_LIB([mat], [matOpen],[],[
dnl    AC_MSG_ERROR(libmat (Mathworks) library is required for this program)])
echo $LDFLAGS

dnl The compilation check for matlabs libs will fail if we dont try to compile
dnl as a shared library.
dnl Also, if third arg to AC_CHECK_LIB is missing, it will automatically append, e.g.,
dnl -lmex to $LIBS, which will make the rest of the checks fail.
tmp_cflags=$CFLAGS
CFLAGS="-shared $CFLAGS"
AC_CHECK_LIB([mex], [mexErrMsgTxt],[
             MATLAB_LIBS=mex],[
   AC_MSG_ERROR(Matlab libmex library not found.)])
CFLAGS=$tmp_cflags

#AC_CHECK_LIB([mwbuiltinsutil], [main])
#AC_CHECK_LIB([mwservices], [main])
#AC_CHECK_LIB([mx], [main])

fi


echo $CFLAGS

AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wextra -Wunused -Werror  \
                         -Wcast-qual -Wcast-align -Winit-self      \
                         -Wfloat-equal -Wundef -Wshadow -Wswitch-enum -Wvla], [WARNFLAGS], [] [])
AC_SUBST(WARNFLAGS, $WARNFLAGS)

if test x$with_opt = xyes;
then
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [OPT], [], [])
fi
AC_SUBST(OPT, $OPT)


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





AM_CONDITIONAL([comp_with_mkl], [test x$use_mkl = xyes])
AM_CONDITIONAL([comp_with_fftw], [test x$use_fftw3 = xyes])



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

AC_MSG_NOTICE(L1C will be built with OPT = $OPT")
AC_MSG_NOTICE(Compile with MKL? $use_mkl)
AC_MSG_NOTICE(Compile with FFTW? $use_fftw3)

AC_CONFIG_FILES([Makefile  src/Makefile src/vcl/Makefile interfaces/Makefile])
AC_OUTPUT
