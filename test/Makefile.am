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
						ABS_TOP_SRCDIR='$(abs_top_srcdir)'; export ABS_TOP_SRCDIR;


TEST_DATA_DIR="$(abs_srcdir)/test_data/"

# -------------------------------------------------------------------------------#
# The valgrind checks are not run by default because they take forever. To do so,
# invoke
# >>> with_valgrind=yes make check
#
$(info "---------------------$(HOSTOS)")
TESTS = test_l1c           \
		build_dct2_data.py     \
		build_dct_data.py     \
		build_cgsolve_data.py \
		build_l1qc_data.py    \
		test_matlab_dlopen.sh \
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
					TEST_l1c_common.c   \
					TEST_l1qc_newton.c  \
					TEST_vcl_math.c     \
					test_l1c.c          \
					TEST_dct.c          \
					TEST_dct2.c

test_l1c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

test_l1c_CPPFLAGS = -D_EXCLUDE_PUBLIC_ -I$(top_srcdir)/include  \
					-I$(top_srcdir)/src \
						 $(L1CCPPFLAGS)
test_l1c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
test_l1c_LDFLAGS = $(LDFLAGS)

log_deps = 	build_dct2_data.log \
		   	build_dct_data.log \
		   	build_cgsolve_data.log \
			build_l1qc_data.log


test_l1c.log:$(log_deps)



test_matlab_dlopen.sh:$(log_deps)
test_matlab_mex.sh:$(log_deps)
test_python_interface.sh:$(log_deps)
test_l1qc_dct_c.sh:$(log_deps)
valgrind_wrapper.log:$(log_deps)
