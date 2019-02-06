#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AM_CFLAGS =	 $(MAT_DEF)

# libl1cinclude_HEADERS =../include/l1qc_newton.h
# Using noinst_LTLLIBRARIES instead of noinst_LIBRARIES prevents the warning
# about linking against a static library when we compile and link the interfaces


noinst_LTLIBRARIES = libl1c.la
libl1c_la_CPPFLAGS = -I$(top_srcdir)/include $(L1CCPPFLAGS)
libl1c_la_CFLAGS = $(L1CCFLAGS)  -fPIC $(WARNFLAGS)
libl1c_la_SOURCES = l1qc_newton.c \
					l1c_math.c \
					cgsolve.c  \
					dct.c      \
					dct_mkl.c


# include_HEADERS = ../include/*.c

# 					../include/cgsolve.h \
# 					../include/l1qc_newton.h    \
# 					../include/l1c_common.h \
# 					../include/l1c_math.h \
# 					../include/dct.h \
# 					../include/dct_mkl.h


AM_CPPFLAGS = $(MKL_INCLUDE)
