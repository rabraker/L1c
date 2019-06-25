#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

# AM_CFLAGS =


lib_LTLIBRARIES = libl1c.la

libl1c_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include $(L1CCPPFLAGS)

libl1c_la_CFLAGS = $(L1CCFLAGS)
libl1c_la_SOURCES = l1qc_newton.c \
					l1c_math.c \
					cgsolve.c  \
					json_utils.c \
					TV.c         \
					bregman.c    \
					l1qc_dct.c   \
					linesearch.c \
					l1c_logging.c \
					l1c_memory.c \
					matrix_transforms.c \
					dct_transforms.c \
                    dct1.c \
					dct2.c \
                    dct_tv_transforms.c \
					nesta.c \
					$(top_srcdir)/include/bregman.h \
					$(top_srcdir)/include/json_utils.h \
					$(top_srcdir)/include/l1c_math.h \
					$(top_srcdir)/include/linesearch.h \
					$(top_srcdir)/include/l1qc_newton.h \
					$(top_srcdir)/include/TV.h \
					$(top_srcdir)/include/l1c_timing.h

libl1c_la_LIBADD	= ../src/vcl/libvcl_math.la

libl1c_la_LDFLAGS	= -no-undefined $(LDFLAGS)\
						$(WIN32LDFLAGS) -Wl,--no-as-needed


