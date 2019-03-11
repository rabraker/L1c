#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies

AM_TESTS_FD_REDIRECT = 9>&1
AM_TESTS_ENVIRONMENT = TOP_BUILDDIR='$(top_builddir)'; export TOP_BUILDDIR; \
						ABS_TOP_BUILDDIR='$(abs_top_builddir)'; export ABS_TOP_BUILDDIR; \
						ABS_TOP_SRCDIR='$(abs_top_srcdir)'; export ABS_TOP_SRCDIR;  \
						WITH_MEX='$(WITH_MEX)'; export WITH_MEX;


TEST_DATA_DIR="$(abs_srcdir)/test_data/"

# -------------------------------------------------------------------------------#
# The valgrind checks are not run by default because they take forever. To do so,
# invoke
# >>> with_valgrind=yes make check
#

TESTS = test_l1c           \
		build_CS20NG_example_data.py \
		build_dct2_data.py     \
		build_dct_data.py     \
		build_cgsolve_data.py \
		build_l1qc_data.py    \
		build_TV_data.py      \
		build_bregman_data.py \
		test_matlab_mex.sh \
		test_python_interface.sh \
		test_l1qc_dct_c.sh       \
		valgrind_wrapper.sh


TEST_EXTENSIONS = .py .sh

# Expected failures
XFAIL_TESTS=

check_PROGRAMS   = test_l1c

test_l1c_SOURCES = TEST_cgsolve.c       \
					TEST_l1c_math.c     \
					TEST_l1c_memory.c   \
					TEST_l1qc_newton.c  \
					TEST_vcl_math.c     \
					TEST_bregman.c      \
					TEST_TV.c           \
					test_l1c.c

if compile_with_mkl
test_l1c_SOURCES +=	TEST_dct.c
else
test_l1c_SOURCES +=	TEST_dct.c TEST_dct2.c
endif

test_l1c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

test_l1c_CPPFLAGS = -D_EXCLUDE_PUBLIC_ -I$(top_srcdir)/include  \
					-I$(top_srcdir)/src \
						 $(L1CCPPFLAGS)
test_l1c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
test_l1c_LDFLAGS = $(LDFLAGS)

log_deps = 	build_CS20NG_example_data.log \
			build_dct2_data.log \
		   	build_dct_data.log \
		   	build_cgsolve_data.log \
			build_l1qc_data.log    \
			build_TV_data.log      \
			build_bregman_data.log


test_l1c.log:$(log_deps)



# test_matlab_dlopen.sh:$(log_deps)

test_matlab_mex.sh:$(log_deps)
test_python_interface.sh:$(log_deps)
test_l1qc_dct_c.sh:$(log_deps)
valgrind_wrapper.log:$(log_deps)
