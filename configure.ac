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
dnl Need this for AC_EGREP_HEADER to work
AC_PROG_EGREP
AC_PROG_LN_S
dnl AC_PROG_MAKE_SET
dnl AC_PROG_CPP
AM_PROG_AR

AC_PROG_INSTALL
LT_INIT([win32-dll, dlopen])
# callint LT_INIT means we dont need: AC_PROG_RANLIB, AM_PROG_LIBTOOL
AC_CANONICAL_HOST
AM_INIT_AUTOMAKE

AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])
dnl This is necessary, otherwise, when we use AC_CHECK_HEADER inside the conditional,
dnl AC_PROG_CPP will be expanded in the else, and not visible outside, which leads to
dnl warnings like WARNING: check.h: accepted by the compiler, rejected by the preprocessor!
AC_PROG_CPP

AH_TEMPLATE([HAVE_SATLAS],
    [Define to 1 if you have libsatlas])
AH_TEMPLATE([HAVE_OPENBLAS],
    [Define to 1 if you have libopenblas])


AC_ARG_WITH(debug, [AS_HELP_STRING([--with-debug],[Compile with debugging support.])],
                   [AC_DEFINE(WITH_DEBUG,1,[Debugging enabled])
                   with_debug=yes], [with_debug=no])

AC_ARG_WITH(opt, [AS_HELP_STRING([--without-opt],[Compile optimizations off.])],
                 [with_opt=no], [with_opt=yes])

AC_ARG_WITH(mkl, [AS_HELP_STRING([--with-mkl], [Compile using Intel Math Kernel Library])], [
                 with_mkl=yes], [with_mkl=no])

AC_ARG_WITH(mex, [AS_HELP_STRING([--with-mex],
                 [Compile matlab interface mex functions. This requires that
                 CPPLAGS contains -I${MATLABROOT}/extern/include and that
                 LDFLAGS contains -L${MATLABROOT}/bin/(glnxa64|win64)/])], [
                 with_mex=yes], [with_mex=no])
AC_ARG_WITH(fftw3, [AS_HELP_STRING([--with-fftw3], [Compile example with FFTW])], [
                 use_fftw3=yes], [use_fftw3=no])
AC_ARG_WITH(fftw3_threads, [AS_HELP_STRING(
                 [--with-fftw3-threads=@<:@combined@:>@@<:@omp@:>@@<:@threads@:>@@<:@yes@:>@@<:@no@:>@],
                   [Compile example with FFTW threading. By default, searches for combined, then omp
                   then threads. The default is to compile with threads. Threading can be
                   disabled by setting ARG=no])], [], [with_fftw3_threads=yes   ])


if test x$with_debug = xyes;
then
    AX_APPEND_COMPILE_FLAGS([-O0 -ggdb -g], [l1c_cflags], [], [])
fi

if test x$with_opt = xyes;
then
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [l1c_cflags], [], [])
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [l1c_cxxflags], [], [])
fi


# --------------------------------------------------------------------------------------- #
# We will build these flags up as we go, and AC_SUBST them at the bottom.


AX_APPEND_ALL_COMPILE_FLAGS([-fopenmp -m64 -std=c11], [l1c_cflags],
                            [AC_MSG_ERROR("Compiler does not support required flag")],
                            [],[])
l1c_cppflags=
mex_cppflags=
mex_cflags=
mkl_cppflags=
mkl_cflags=
test_cppflags=
test_cflags=
win32_ldflags=

# At least on my linux, this flag is not supported by ld, but is needed to
# build on mingw
AX_APPEND_LINK_FLAGS([-Wl,--export-all], [win32_ldflags], [], [])

AC_CHECK_LIB([m], [sqrt], [], [])


# Check for posix_memalign or _mm_malloc
AC_CHECK_FUNC(posix_memalign, [have_posix_mem=yes], [])

AC_LANG([C])

# ------------------------------ How can we align memory?----------------------------------- #
# _mm_malloc is an inline function in mm_malloc.h, included by xmmintrin.h
# Should not need a header check, since if xmmintrin.h is not found, this should fail.
AC_MSG_CHECKING([checking for _mm_malloc])
AC_LINK_IFELSE([AC_LANG_SOURCE([[
#include <xmmintrin.h>
void main(){
     double *x=_mm_malloc(10*sizeof(double), 64);
}]])], [have_mm_malloc=yes
       AC_MSG_RESULT([yes])], [AC_MSG_RESULT([no])])

if test x$have_posix_mem = xyes;
then
    AC_DEFINE([_POSIX_C_SOURCE], [200112L], [Unlocks posix_memalign in stdlib.h.])
    AC_DEFINE([_HAVE_POSIX_MEMALIGN_], [1], [We have posix_memalign.])

elif test x$have_mm_malloc = xyes;
then
     AC_DEFINE(_HAVE_MM_MALLOC_, [1], [We have _mm_malloc])
else
    AC_MSG_ERROR(This program requires either posix_memalign() or _mm_malloc in the standard library.);
fi


# ----------------------------------------------------------------------------------------- #
# ------------------------------ Checks for MKL ------------------------------------------- #
if test "x$with_mkl" = "xyes";
then
# Intel says the link order should look like:
# -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
# AC_CHECK_LIB prepends the libraries, so we work backwards.
# N.B. If we dont also linke with libkl_avx, we get link errors when dlopening in python and matlab.
# Intel has nothing to say about this in the link advisor...
    AC_CHECK_LIB([dl], [dlopen],[], [AC_MSG_ERROR("Cannot find libdl")])

    AC_CHECK_LIB([pthread], [pthread_mutex_trylock], [],
                 [AC_MSG_ERROR([Cannot find libpthread, which is required
                 with --with-fftw3-threads=threads])])

    AC_CHECK_LIB([gomp], [omp_set_num_threads],[],[
                 AC_MSG_ERROR(omp library is required for this program)])

    mkl_libs_check="-lmkl_core -lmkl_gnu_thread -lmkl_avx -lgomp"
    AC_SEARCH_LIBS([mkl_xblas_avx_blas_ztrmv_d], [mkl_avx], [],
               [AC_MSG_ERROR("Cannot find libmkl_core")], [$mkl_libs_check])
    AC_SEARCH_LIBS([mkl_spblas_zzeros], [mkl_core], [],
               [AC_MSG_ERROR("Cannot find libmkl_core")], [$mkl_libs_check])
    AC_SEARCH_LIBS([mkl_trs_djacobi], [mkl_gnu_thread], [],
               [AC_MSG_ERROR("Cannot find libmkl_gnu_thread")], [$mkl_libs_check])
    AC_SEARCH_LIBS([cblas_daxpby], [mkl_intel_lp64], [],
               [AC_MSG_ERROR("Cannot find libmkl_intel_lp64")], [$mkl_libs_check])

    AC_CHECK_HEADER([mkl.h], [], [AC_MSG_ERROR("Cannot find mkl.h")])

    AC_DEFINE([_USEMKL_], [1], [We have MKL.])
    AC_DEFINE([MKL_LP64], [1], [Use ILP64 with mkl.])

else
 dnl ----------- look for another BLAS library ---------------------
 AC_CHECK_HEADER([cblas.h], [], [AC_MSG_ERROR("Cannot find cblas.h")])

    if test ! -z $BLAS_LIB ; then
       AC_CHECK_LIB([$BLAS_LIB], [cblas_daxpy],[
                    tmp="HAVE_${BLAS_LIB^^}"
                    AC_DEFINE_UNQUOTED([$tmp], [1])
                    LIBS="-l$BLAS_LIB $LIBS"], [BLAS_LIB=])
    fi
    # Check in OpenBlas
    if test -z $BLAS_LIB; then
       AC_CHECK_LIB([openblas], [cblas_daxpy],[BLAS_LIB=openblas
                    AC_DEFINE([HAVE_OPENBLAS], [1])
                    LIBS="-lopenblas $LIBS"],[])
    fi
    # Check in ATLAS
    if test -z $BLAS_LIB; then
       AC_CHECK_LIB([satlas], [cblas_daxpy],[BLAS_LIB=satlas
       AC_DEFINE([HAVE_SATLAS], [1], [0])
       LIBS="-lsatlas $LIBS"],[])
    fi

    # If we still havent found it, error out.
    if test -z $BLAS_LIB; then
       AC_MSG_ERROR(BLAS library with cblas_daxpy is required for this program)
    fi

# Now check if the supplied lib has daxpby
  case $BLAS_LIB in
       openblas)
            AC_CHECK_FUNC([cblas_daxpby], [AC_DEFINE([HAVE_DAXPBY], [1], [0])], [have_daxpby=no])
            AC_MSG_CHECKING(checking for cblas_daxpby in cblas.h)
            AC_EGREP_HEADER([cblas_daxpby], [cblas.h], [AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)
                          AC_MSG_ERROR(Could not find cblas_daxpby in cblas.h)])
            ;;
       satlas)  # Make a special exception for ATLAS, since they call it catlas_daxpby
            AC_CHECK_FUNC([catlas_daxpby], [
            AC_DEFINE([HAVE_DAXPBY], [1], [0])],
            [have_daxpby=no])
            AC_MSG_CHECKING(checking for cblas_daxpby in cblas.h)
            AC_EGREP_HEADER([catlas_daxpby], [cblas.h], [AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)
                           AC_MSG_ERROR(Could not find cblas_daxpby in cblas.h)])
            ;;
           *)
           AC_CHECK_FUNC(cblas_daxpby, [AC_DEFINE([HAVE_DAXPBY], [1], [0])], [have_daxpby=no])
           AC_MSG_CHECKING(checking for cblas_daxpby in cblas.h)
           AC_EGREP_HEADER([cblas_daxpby], [cblas.h], [AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)
                           AC_MSG_ERROR(Could not find cblas_daxpby in cblas.h)])
           ;;
  esac

  #For now, throw an error if we dont have cblas_daxpby. Later, we will provide an implementation.
  if test x$have_daxpby = xno; then
     AC_MSG_ERROR(BLAS library with cblas_daxpby is required for this program)
  fi

fi

# ------------------------------ Checks for FFTW ---------------------------------------- #
if test x$use_fftw3 = xyes;
then
   AC_CHECK_HEADER([fftw3.h], [], [AC_MSG_ERROR("Cannot find fftw3.h")])

   AC_CHECK_LIB([fftw3], [fftw_plan_r2r_1d],[],[
                         AC_MSG_ERROR(fftw3 library is required for this program)])
# Figure out which threading library we are supposed to use, if any.
   have_fftw_threads=no
   maybe_suffx="combined omp threads"
   for suffx in $maybe_suffx; do
       case $suffx in
            combined)
                lib=fftw3
                ;;
            omp)
                lib=fftw3_omp
                ;;
            threads)
                lib=fftw3_threads
                ;;
            *)
       esac

       if test x$suffx = x$with_fftw3_threads ||test x$with_fftw3_threads = xyes   ; then
          AC_CHECK_LIB([$lib], [fftw_init_threads],[have_fftw_threads=yes],[have_fftw_threads=no])
       fi

       if test x$have_fftw_threads = xyes; then
          LIBS="$LIBS -l$lib"
          break
       fi
   done

   if test x$with_fftw3_threads = xyes && test x$have_fftw_threads = xno; then
      AC_MSG_ERROR(fftw3 threading library not found)]
   elif test x$have_fftw_threads = xyes; then
      AC_DEFINE([HAVE_FFTW3_THREADS], [1], [Compile with fftw threading support])
   fi

   # We link with OMP anyway, down below. But if the user selected thread,
   # we need to also link with pthread.
   if test x$lib = xthread; then
   AC_CHECK_LIB([pthread], [pthread_mutex_trylock], [],
                 AC_MSG_ERROR([Cannot find libpthread, which is required with --with-fftw3-threads=threads]))
   fi

  AC_DEFINE([_USEFFTW3_], [1], 0)
fi  # checks for fftw

echo "Compile with fftw3? $use_fftw3"
echo "Compile with fftw3 threads? $have_fftw_threads"



dnl -Werror
AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wextra -Wunused   \
                          -Wcast-qual -Wcast-align -Winit-self     \
                          -Wfloat-equal -Wundef -Wshadow -Wstringop-truncation \
                          -Wswitch-enum -Wvla], [WARNFLAGS], [] [])
AC_SUBST(WARNFLAGS, $WARNFLAGS)

# ----------------------------------------------------------------------------------------- #
# ------------------------------ Checks for matlab MEX------------------------------------- #

if test x$with_mex = xyes;
then

dnl The compilation check for matlabs libs will fail if we dont try to compile
dnl as a shared library.
dnl Also, if third arg to AC_CHECK_LIB is missing, it will automatically append, e.g.,
dnl -lmex to $LIBS, which will make the rest of the checks fail.
    tmp_cflags=$CFLAGS
    CFLAGS="-shared $CFLAGS"
    AC_CHECK_LIB([mex], [mexErrMsgTxt],[mex_ldlibs="-lmex"],[
                 AC_MSG_ERROR(Matlab libmex library not found.)])
    AC_CHECK_LIB([mx], [mxIsScalar], [mex_ldlibs="$mex_ldlibs -lmx"],
                 [AC_MSG_ERROR(Matlab libmx library not found.)])
    AC_CHECK_LIB([mwservices], [svMatlabRoot], [mex_ldlibs="$mex_ldlibs -lmwservices"],
                 [AC_MSG_ERROR(Matlab libmwservices library not found.)])
    # This lib is not present on windows but seems to build fine. Should check
    # if it is really needed on linux.
    AC_CHECK_LIB([mwbuiltinsutil], [_fini], [mex_ldlibs="$mex_ldlibs -lmwbuiltinsutil"],
                 [])
# Reset the user CFLAGS
    CFLAGS=$tmp_cflags

    AC_DEFINE([__MATLAB__], [1], [1 to compile with mex interfaces])
              mex_cflags="-fexceptions"

    if test -z $MEXEXT; then
       case $host_os in
            mingw*) MEXEXT=$(mexext.bat) ;;
            *) MEXEXT=$(mexext) ;;
       esac
       if test $? -ne 0; then
          AC_MSG_ERROR([Mathworks mexext must be on the system path, or you must export
                    the proper mex file extension as the envrionmental variable MEXEXT.])
       fi
    fi
    AC_SUBST(MEXEXT, [.$MEXEXT])
fi  # checks for mex


# Checks for programs.

# Checks for libraries.
AC_CHECK_LIB([check], [suite_add_tcase],[],[
             AC_MSG_ERROR(check library is required for this program)])
AC_CHECK_LIB([cjson], [cJSON_GetArraySize],[],[
             AC_MSG_ERROR(cJSON library is required for this program)])

AC_CHECK_LIB([gomp], [omp_set_num_threads],[],[
             AC_MSG_ERROR(omp library is required for this program)])

AC_CHECK_LIB([m], [sqrt],[],[
             AC_MSG_ERROR(library is required for this program)])


# Check for headers

AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])
AC_CHECK_HEADER([check.h], [], [AC_MSG_ERROR("Cannot find check.h")], [])
AC_CHECK_HEADER([cjson/cJSON.h], [], [AC_MSG_ERROR("Cannot find CJSON.h")])
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR("Cannot find math.h")])

AM_CONDITIONAL([compile_dct_example], [test x$with_mkl = xyes || test x$use_fftw3 = xyes])
AM_CONDITIONAL([compile_with_mkl], [test x$with_mkl = xyes ])
AM_CONDITIONAL([comp_with_mex], [test x$with_mex = xyes])


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

# ----------------------------------------------------------------- #
# Make the flags we collected visible to automake
AC_SUBST(L1CCXXFLAGS, $l1c_cxxflags)

AC_SUBST(L1CCPPFLAGS, $l1c_cppflags)
AC_SUBST(L1CCFLAGS, $l1c_cflags)

AC_SUBST(MEXCPPFLAGS, $mex_cppflags)
AC_SUBST(MEXCFLAGS, $mex_cflags)
AC_SUBST(MEXLDLIBS, $mex_ldlibs)
AC_SUBST(WITH_MEX, $with_mex)
AC_SUBST(MKLCPPFLAGS, $mkl_cppflags)
AC_SUBST(MKLCFLAGS, $mkl_cflags)

AC_SUBST(TESTCPPFLAGS, $test_cppflags)
AC_SUBST(TESTCFLAGS, $test_cflags)
AC_SUBST(WIN32LDFLAGS, $win32_ldflags)

AC_MSG_NOTICE(L1C will be built with l1c_cflags = $l1c_cflags")
AC_MSG_NOTICE(L1C will be built with l1c_cppflags = $l1c_cppflags")
AC_MSG_NOTICE(Compile with MKL? $with_mkl)
AC_MSG_NOTICE(Compile with FFTW? $use_fftw3)
AC_MSG_NOTICE(Using mex linkage $mex_ldlibs)

AC_CONFIG_FILES([Makefile  test/Makefile src/Makefile src/vcl/Makefile interfaces/Makefile])
AC_OUTPUT
