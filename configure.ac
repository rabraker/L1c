#                                               -*- Autoconf -*-
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
AM_INIT_AUTOMAKE([-Wall])

AM_CONFIG_HEADER([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Enable using pkg-config. First check if pkg-config is installed.
AC_CHECK_PROG(PKG_CONFIG,[pkg-config],[pkg-config],[no])
test "$PKG_CONFIG" == "no" && AC_MSG_ERROR([Required program 'pkg-config' not found.])

PKG_PROG_PKG_CONFIG

dnl This is necessary, otherwise, when we use AC_CHECK_HEADER inside the conditional,
dnl AC_PROG_CPP will be expanded in the else, and not visible outside, which leads to
dnl warnings like WARNING: check.h: accepted by the compiler, rejected by the preprocessor!
AC_PROG_CPP
AC_LANG([C])

AC_ARG_VAR(BLAS_LIB, [BLAS_LIB=[openblas][satlas][tatlas][?], the name of your prefered BLAS library.
                     By default, configure will search first for openblas, then tatlas, then satlas])

AC_ARG_VAR(python_prefix, [python_prefix=/path/to/python-binding/installation,
                          The path to which the python bindings should be installed.
                          This will default to the path provided by distutils.
                          This variable has no effect if python bindings are not being built.])

AC_ARG_VAR(mex_prefix, [mex_prefix=/path/to/mex/installation,
                       Location to install the mex bindings in. Ignored if --with-mex is
                       not also given. By default, the bindings will be installed in ${prefix}/lib.])

# Needed for AC_DEFINE_UNQUOTED below.
dnl AH_TEMPLATE([HAVE_ATLAS],
dnl     [Define to 1 if you have libsatlas])
dnl AH_TEMPLATE([HAVE_OPENBLAS],
dnl     [Define to 1 if you have libopenblas])

AC_ARG_WITH(debug, [AS_HELP_STRING([--with-debug],[Compile with debugging support.])],
                   [AC_DEFINE(WITH_DEBUG,1,[Debugging enabled])
                   with_debug=yes], [with_debug=no])

AC_ARG_WITH(opt, [AS_HELP_STRING([--without-opt],[Compile optimizations off.])],
                 [with_opt=no], [with_opt=yes])

AC_ARG_WITH(fftw3_threads, [AS_HELP_STRING(
                 [--with-fftw3-threads=@<:@combined@:>@@<:@omp@:>@@<:@threads@:>@@<:@yes@:>@@<:@no@:>@],
                   [Compile example with FFTW threading. By default, searches for combined, then omp
                   then threads. The default is to compile with threads. Threading can be
                   disabled by setting ARG=no])], [], [with_fftw3_threads=yes   ])

AC_ARG_ENABLE(mex, [AS_HELP_STRING([--enable-mex],
                 [Compile matlab interface mex bindings. This requires that
                 CPPLAGS contains -I${MATLABROOT}/extern/include and that
                 LDFLAGS contains -L${MATLABROOT}/bin/(glnxa64|win64)/])], [
                 enable_mex=yes], [enable_mex=no])

AC_ARG_ENABLE([python],[AS_HELP_STRING([--enable-python],
                [build with python bindings])],
                [enable_python=$enableval],[enable_python=no])

if test x$mex_prefix = x; then
   mex_prefix=$libdir
else
# Remove any trailing slash the user may have provided,
# because this breaks the install with libtool, and will result in, e.g.,
#  libtool: error: cannot install 'l1qc_dct_mex.la' to a directory not ending in /home/arnold/matlab/l1c-mex-test/.
    mex_prefix=`echo $mex_prefix|sed 's,\/$,,'`
fi

# --------------------------------------------------------------------------------------- #
# We will build these flags up as we go, and AC_SUBST them at the bottom.
l1c_test_libs=
l1c_cppflags=
l1c_cflags=
l1c_cxxflags=
mex_cppflags=
mex_cflags=
test_cppflags=
test_cflags=
win32_ldflags=

if test x$with_debug = xyes;
then
    AX_APPEND_COMPILE_FLAGS([-O0 -ggdb -g], [l1c_cflags], [], [])
    AC_LANG([C++])
    AX_APPEND_COMPILE_FLAGS([-O0 -ggdb -g], [l1c_cxxflags], [], [])
    AC_LANG([C])
fi

if test x$with_opt = xyes;
then
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [l1c_cflags], [], [])
    AC_LANG([C++])
    AX_APPEND_COMPILE_FLAGS([-O3 -march=native -mtune=native], [l1c_cxxflags], [], [])
    AC_LANG([C])
fi



AX_APPEND_ALL_COMPILE_FLAGS([-fopenmp -m64 -std=c11], [l1c_cflags],
                            [AC_MSG_ERROR("Compiler does not support required flag")],
                            [],[])

# Check if c compiler accepts a bunch of warning flags.
# Force it to fail on a warning.
AX_CHECK_COMPILE_FLAG(-Werror, [tmp_cflags="-Werror"], [], [], [])
AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wextra -Wunused   \
                            -Wcast-qual -Wcast-align -Winit-self     \
                            -Wfloat-equal -Wundef -Wshadow -Wstringop-truncation \
                            -Wswitch-enum -Wvla], [l1c_cflags], [$tmp_cflags] [])

AC_LANG([C++])
AX_CHECK_COMPILE_FLAG(-Werror, [tmp_cflags="-Werror"], [], [], [])
AX_APPEND_ALL_COMPILE_FLAGS([-m64 -std=c++11], [l1c_cxxflags],
                            [AC_MSG_ERROR("Compiler does not support required flag")],
                            [],[])
# Check if c compiler accepts a bunch of warning flags
AX_APPEND_COMPILE_FLAGS([-pedantic -Wall -Wextra -Wunused   \
                            -Wcast-qual -Wcast-align -Winit-self     \
                            -Wfloat-equal -Wundef -Wshadow -Wstringop-truncation \
                            -Wswitch-enum -Wvla], [l1c_cxxflags], [$tmp_cflags] [])
AC_LANG([C])


# At least on my linux, this flag is not supported by ld, but is needed to
# build on mingw
AX_APPEND_LINK_FLAGS([-Wl,--export-all], [win32_ldflags], [], [])



# ----------------------------------------------------------------------------------------- #
# ---------------------------- Checks building python bindings ---------------------------- #
VALID_PYTHON=
if test x$enable_python != xno; then
   if test ! -z $python_prefix; then
      # These variables only have an affect on AX_PYTHON_DEVEL, no
      # AM_PATH_PYTHON.
      PYTHON_SITE_PKG_EXEC=$python_prefix
      PYTHON_SITE_PKG=$python_prefix
   fi

   AM_PATH_PYTHON([3], [PYTHON_FOUND=true],[])
   if test $PYTHON_FOUND = true; then
      AX_PYTHON_DEVEL([>= '3'])
      if test "x$PYTHON" = x ; then
         AC_MSG_ERROR([Incomplete python development package])
         VALID_PYTHON=no
      else
        VALID_PYTHON=yes
      fi

      # We want to install the bindings and their related python files all in the
      # same place. Automake will try to stick these in different directories. And the
      # default will be /usr/lib/python3.5/...
      pyexecdir="${PYTHON_SITE_PKG_EXEC}/${PACKAGE}Py"
      pkgpythondir=$pyexecdir
      pythondir=$pyexecdir
      if test -z $pyexecdir;
      then
        AC_MSG_ERROR(Could not locate PYTHON_EXEC_DIR)
      fi

# look for the numpy include directory. DO NOT INDENT THIS.
npy_cppflags=`$PYTHON -c "import numpy.distutils as npd; \
inc_dirs='';
for dir in npd.misc_util.get_numpy_include_dirs():
    inc_dirs = ' -I'+dir;
print(inc_dirs);"`

      # Other Variables that may prove useful:
      # $pyexecdir, $pythondir, $PYTHON_SITE_PKG, $PYTHON_SITE_PKG_EXEC
      # $npy_cppflags, $enable_python, $PYTHON, $PYTHON_VERSION
      # $PYTHON_CPPFLAGS, $PYTHON_LIBS, $PYTHON_EXTRA_LDFLAGS
      dnl AC_SUBST(PYTHONDIR,$pythondir)
      dnl AC_SUBST(PKGPYTHONDIR, $pkgpythondir)
      dnl AC_SUBST(PYEXECDIR, $pyexecdir)
      dnl AC_SUBST(PKGPYEXECDIR, $pkgpyexecdir)
      AC_SUBST(NPYCPPFLAGS, $npy_cppflags)
   fi
fi


AM_CONDITIONAL([ENABLE_PYTHON],[test x$VALID_PYTHON = xyes])


AC_CHECK_LIB([m], [sqrt], [], [])

# ------------------------------ How can we align memory?----------------------------------- #
# Check for posix_memalign or _mm_malloc
AC_CHECK_FUNC(posix_memalign, [have_posix_mem=yes], [])


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



dnl ----------- look for a BLAS library ---------------------
AC_CHECK_HEADER([cblas.h], [], [AC_MSG_ERROR("Cannot find cblas.h")])
# All blas libraries should cblas_daxpy. Not all have cblas_daxpby. Check
# cblas_daxpy first.
if test ! -z $BLAS_LIB; then
   AC_CHECK_LIB([$BLAS_LIB], [cblas_daxpy],[
                tmp="HAVE_${BLAS_LIB^^}"
                LIBS="-l$BLAS_LIB $LIBS"], [BLAS_LIB=])
fi
# Check in OpenBlas
if test -z $BLAS_LIB; then
   AC_CHECK_LIB([openblas], [cblas_daxpy],[BLAS_LIB=openblas
                LIBS="-lopenblas $LIBS"],[])
fi
# Check in ATLAS
if test -z $BLAS_LIB; then
# Atlas can compile both a serial (libsatlas) and threaded (libtatlas), library.
# Check the threaded library first.
   AC_CHECK_LIB([tatlas], [cblas_daxpy],[BLAS_LIB=tatlas
                LIBS="-ltatlas $LIBS"],[
                          AC_CHECK_LIB([satlas], [cblas_daxpy],[BLAS_LIB=satlas
                          LIBS="-lsatlas $LIBS"],[])])
fi
# Check in system blas
if test -z $BLAS_LIB; then
AC_CHECK_LIB([blas], [cblas_daxpy],[BLAS_LIB=blas
             LIBS="-lblas $LIBS"],[])
fi
# If we still havent found it, error out.
if test -z $BLAS_LIB; then
   AC_MSG_ERROR(BLAS library with cblas_daxpy is required for this program)
fi

# Now check if the supplied lib has daxpby
case $BLAS_LIB in
  satlas|tatlas)  # Make a special exception for ATLAS, since they call it catlas_daxpby
       AC_CHECK_FUNC([catlas_daxpby], [have_daxpby=yes], [have_daxpby=no])
       if test x$have_daxpby = xyes; then
          AC_MSG_CHECKING(checking for catlas_daxpby in cblas.h)
          AC_EGREP_HEADER([catlas_daxpby], [cblas.h], [
          AC_MSG_RESULT(yes)
          AC_DEFINE([HAVE_ATLAS_DAXPBY], [1], [0])], [AC_MSG_RESULT(no)])
       fi
       ;;
       *)
       AC_CHECK_FUNC([cblas_daxpby], [have_daxpby=yes], [have_daxpby=no])
       if test x$have_daxpby = xyes; then
          AC_MSG_CHECKING(checking for cblas_daxpby in cblas.h)
          AC_EGREP_HEADER([cblas_daxpby], [cblas.h], [
          AC_MSG_RESULT(yes)
          AC_DEFINE([HAVE_CBLAS_DAXPBY], [1], [0])], [AC_MSG_RESULT(no)])
       fi
       ;;
esac

dnl fi

# ------------------------------ Checks for FFTW ---------------------------------------- #
AC_CHECK_HEADER([fftw3.h], [], [AC_MSG_ERROR("Cannot find fftw3.h")])

AC_CHECK_LIB([fftw3], [fftw_plan_r2r_1d],[],[
                      AC_MSG_ERROR(fftw3 library is required for this program)])
_fftw_ldlibs="-lfftw3"
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
    # If user specified a specific threading library,
    # skip until the loop reaches it.
    if test x$suffx = x$with_fftw3_threads ||test x$with_fftw3_threads = xyes   ; then
       AC_CHECK_LIB([$lib], [fftw_init_threads],[have_fftw_threads=yes],[have_fftw_threads=no])
    fi

    if test x$have_fftw_threads = xyes; then
       LIBS="$LIBS -l$lib"
       _fftw_ldlibs="-l$lib $_fftw_ldlibs"
       break
    fi
done

if test ! x$with_fftw3_threads = xno && test x$have_fftw_threads = xno; then
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

echo "Compile with fftw3 threads? $have_fftw_threads"




# ----------------------------------------------------------------------------------------- #
# ------------------------------ Checks for matlab MEX------------------------------------- #

if test x$enable_mex = xyes;
then

dnl The compilation check for matlabs libs will fail if we dont try to compile
dnl as a shared library.
dnl Also, if third arg to AC_CHECK_LIB is missing, it will automatically append, e.g.,
dnl -lmex to $LIBS, which will make the rest of the checks fail.
    tmp_cflags=$CFLAGS
    CFLAGS="-shared -fPIC $CFLAGS"
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

case $host_os in
     *linux*) SHRLIB_EXT="so"
     ;;
     *darwin*) SHRLIB_EXT="dylib"
     ;;
     *mingw*) SHRLIB_EXT="dll"
     ;;
     *) SHRLIB_EXT="so"
     ;;
esac
echo "====================================="
echo $host_os
echo $host
echo $SHRLIB_EXT
AC_SUBST(SHRLIB_EXT, $SHRLIB_EXT)

# ---------------------- Check for headers -------------------------------

AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h])
AC_CHECK_HEADER([omp.h], [], [AC_MSG_ERROR("Cannot find omp.h")], [])
AC_CHECK_HEADER([check.h], [], [AC_MSG_ERROR("Cannot find check.h")], [])
# We need to make sure the installed libcheck we found has the floating point macros,
# which were introduced in check ver 0.11.0
AC_MSG_CHECKING(checking for ck_assert_double_eq_tol in check.h)
# Cant use AC_EGREP_HEADER, because it greps the output of pre-processor on check.h, not check.h itself.
#
AC_COMPILE_IFELSE([AC_LANG_SOURCE([[
#include "check.h"
#ifndef ck_assert_double_eq_tol
# error Macro ck_assert_double_eq_tol not defined
#endif
]])], [AC_MSG_RESULT(yes)],
      [AC_MSG_RESULT(no)
      AC_DEFINE([HAVE_OLD_LIBCHECK], 1, [Define to 1 if libcheck does not include floating point test macros.])])

AC_CHECK_HEADER([cjson/cJSON.h], [], [AC_MSG_ERROR("Cannot find CJSON.h")])
AC_CHECK_HEADER([math.h], [], [AC_MSG_ERROR("Cannot find math.h")])



# Checks for programs.

# ----------------------- Checks for libraries.---------------------------------
AC_CHECK_LIB([m], [sqrt],[],[
             AC_MSG_ERROR(Math library is required for this program)])


# If libcheck is installed from source, dont need to link with all these extra libs, so we shouldnt require it. Easiest to just use pkg-config.

PKG_CHECK_MODULES([check], [check])

AC_SUBST(check_CFLAGS, $check_CFLAGS)
AC_SUBST(check_LIBS, $check_LIBS)


AC_CHECK_LIB([cjson], [cJSON_GetArraySize],[],[
             AC_MSG_ERROR(cJSON library is required for this program)])

# Try both libgomp and libomp. When compiling with clang, this will artificially pass
# w/o -fopenmp. When we compile L1c, -fopenmp will cause clang to add -lomp to the ld
# call, which will fail, even if we found gomp here.
AC_CHECK_LIB([gomp], [omp_set_num_threads],[],
                     [AC_CHECK_LIB([omp], [omp_set_num_threads],[],
                     [AC_MSG_ERROR(omp library is required for this program)], [-fopenmp])], [-fopenmp])


AM_CONDITIONAL([comp_with_mex], [test x$enable_mex = xyes])


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
AC_SUBST(L1CCFLAGS, $l1c_cflags)


AC_SUBST(MEXCPPFLAGS, $mex_cppflags)
AC_SUBST(MEXCFLAGS, $mex_cflags)
AC_SUBST(MEXLDLIBS, $mex_ldlibs)
AC_SUBST(WITH_MEX, $enable_mex)
AC_SUBST(MEX_PREFIX, $mex_prefix)

AC_SUBST(WITH_PYTHON, $enable_python)


AC_SUBST(TESTCPPFLAGS, $test_cppflags)
AC_SUBST(TESTCFLAGS, $test_cflags)
AC_SUBST(L1CTESTLIBS, $l1c_test_libs)

AC_SUBST(WIN32LDFLAGS, $win32_ldflags)


AC_CONFIG_FILES([Makefile
                test/Makefile
                src/Makefile
                src/vcl/Makefile
                examples/c/Makefile
                include/Makefile
                include/vcl/Makefile
                interfaces/python/Makefile
                interfaces/mex/Makefile
                interfaces/mex/test/Makefile])



AC_OUTPUT([test/init_matlab_paths.m])

echo "  "
echo "================================================================================================="
AC_MSG_NOTICE(Configuration summary:)
AC_MSG_NOTICE(L1C will be built with l1c_cflags = $l1c_cflags)
AC_MSG_NOTICE(L1C will be built with l1c_cxxflags = $l1c_cxxflags)
AC_MSG_NOTICE(L1C will be built with l1c_cppflags = $l1c_cppflags)
AC_MSG_NOTICE(Compile with FFTW threads? $have_fftw_threads: $_fftw_ldlibs)
AC_MSG_NOTICE(Compile python bindings? $enable_python)
AC_MSG_NOTICE(Compile mex bindings? $enable_mex)

AC_MSG_NOTICE(libl1c will be installed into: ${prefix})
if test x$enable_python != xno; then
   AC_MSG_NOTICE(Python bindings will be installed in: $pyexecdir)
   AC_MSG_NOTICE(Python bindings will be installed in: $pkgpythondir)
fi
if test x$enable_mex != xno; then
   AC_MSG_NOTICE(Mex bindings will be installed in: $mex_prefix)
   AC_MSG_NOTICE(Using mex linkage $mex_ldlibs)
fi
