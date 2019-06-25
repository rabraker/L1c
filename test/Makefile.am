#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8


test_data_dir = $(abs_builddir)/test_data

AM_TESTS_FD_REDIRECT = 9>&1
AM_TESTS_ENVIRONMENT = 	ABS_TOP_SRCDIR='$(abs_top_srcdir)'; export ABS_TOP_SRCDIR;  \
						ABS_TOP_BUILDDIR='$(abs_top_builddir)'; export ABS_TOP_BUILDDIR; \
						TEST_DATA_DIR="$(test_data_dir)"; export TEST_DATA_DIR; \
						WITH_MEX='$(WITH_MEX)'; export WITH_MEX;                    \
						WITH_PYTHON='$(WITH_PYTHON)'; export WITH_PYTHON;


# -------------------------------------------------------------------------------#
# The valgrind checks are not run by default because they take forever. To do so,
# invoke
# >>> with_valgrind=yes make check
#

TESTS = test_l1c.test                  \
		build_CS20NG_example_data.py   \
		build_dct2_data.py             \
		build_dct_data.py              \
		build_dctTV_data.py  		   \
		build_cgsolve_data.py          \
		build_l1qc_data.py             \
		build_nesta_data.py            \
		build_TV_data.py               \
		build_bregman_data.py          \
		build_matrix_transform_data.py \
		test_matlab_mex.sh             \
		test_python_interface.sh       \
		test_l1qc_dct_c.sh             \
		valgrind_wrapper.sh


TEST_EXTENSIONS = .py .sh .test

EXTRA_DIST =	build_CS20NG_example_data.py   \
				build_dct2_data.py 	   		   \
				build_dct_data.py  	  		   \
				build_dctTV_data.py  		   \
				build_cgsolve_data.py 		   \
				build_l1qc_data.py 	  		   \
				build_nesta_data.py            \
				build_TV_data.py   	  		   \
				build_bregman_data.py 		   \
				build_matrix_transform_data.py \
				test_matlab_mex.sh 			   \
				test_python_interface.sh 	   \
				test_l1qc_dct_c.sh 		 	   \
				valgrind_wrapper.sh			   \
				test_runner_utils.sh           \
				L1cTestDataUtils.py            \
				check_gradient.py


# Expected failures
XFAIL_TESTS=

check_PROGRAMS   = test_l1c.test

test_l1c_test_SOURCES = TEST_cgsolve.c       \
						TEST_l1c_math.c     \
						TEST_l1c_memory.c   \
						TEST_l1qc_newton.c  \
						TEST_vcl_math.c     \
						TEST_bregman.c      \
						TEST_TV.c           \
						TEST_linesearch.c        \
						TEST_nesta.c             \
						TEST_matrix_transforms.c \
						test_l1c.c               \
						TEST_dct.c               \
						TEST_dct2.c              \
						TEST_dct_tv.c            \
						$(top_srcdir)/include/bregman.h        \
						$(top_srcdir)/include/check_utils.h    \
						$(top_srcdir)/include/json_utils.h     \
						$(top_srcdir)/include/l1c.h            \
						$(top_srcdir)/include/l1c_math.h       \
						$(top_srcdir)/include/l1c_timing.h     \
						$(top_srcdir)/include/l1qc_newton.h    \
						$(top_srcdir)/include/linesearch.h     \
						$(top_srcdir)/include/nesta.h          \
						$(top_srcdir)/include/test_constants.h \
						$(top_srcdir)/include/TEST.h           \
						$(top_srcdir)/include/TV.h             \
						$(top_srcdir)/include/vcl_math.h



test_l1c_test_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la $(check_LIBS)
test_l1c_test_CPPFLAGS = -D_EXCLUDE_PUBLIC_ -I$(top_srcdir)/include  \
						-I$(top_srcdir)/src $(L1CCPPFLAGS)

test_l1c_test_CFLAGS = $(L1C_TEST_CFLAGS) $(check_CFLAGS)
test_l1c_test_LDFLAGS = $(LDFLAGS)

log_deps = 	build_CS20NG_example_data.log \
			build_dct2_data.log \
		   	build_dct_data.log \
			build_dctTV_data.log   \
		   	build_cgsolve_data.log \
			build_nesta_data.log   \
			build_l1qc_data.log    \
			build_TV_data.log      \
			build_bregman_data.log \
			build_matrix_transform_data.log



$(EXTRA_DIST):$(test_data_dir)

$(test_data_dir):
	mkdir -p $(test_data_dir)

# Ensure data gets built first
test_l1c.log:$(log_deps)

test_matlab_mex.log:$(log_deps)

test_python_interface.log:$(log_deps)

test_l1qc_dct_c.log:$(log_deps)

valgrind_wrapper.log:$(log_deps)

distclean-local:
	rm -f $(test_data_dir)/*.json
