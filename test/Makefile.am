#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies



$(info "MATHLIBS = $(MATHLIBS)")
$(info "LDFLAGS = $(LDFLAGS)")
$(info "MEXCPP = $(MEXCPPFLAGS)")
$(info "MEXLIBS = $(MEXLDLIBS)")
# ----------------------------------------------------------------------#
# -------------------------------  ----------------------------------#
# if compile_test

bin_PROGRAMS   = test_l1c

test_l1c_SOURCES = TEST_cgsolve.c      \
					TEST_l1c_math.c     \
					TEST_l1c_common.c   \
					TEST_l1qc_newton.c  \
					TEST_vcl_math.c     \
					test_l1magic.c      \
					TEST_dct.c          \
					TEST_dct_mkl.c
test_l1c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

test_l1c_CPPFLAGS =-I$(top_srcdir)/include  \
						 $(L1CCPPFLAGS)
test_l1c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
test_l1c_LDFLAGS = $(LDFLAGS)


# endif
