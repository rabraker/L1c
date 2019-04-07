#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies

bin_PROGRAMS =l1qc_dct_c
l1qc_dct_c_SOURCES = l1qc_dct_c.c
l1qc_dct_c_LDADD =  $(top_builddir)/src/libl1c.la \
					$(top_builddir)/src/vcl/libvcl_math.la

l1qc_dct_c_CPPFLAGS =-I$(top_srcdir)/include  \
						 $(L1CCPPFLAGS)
l1qc_dct_c_CFLAGS = $(L1CCFLAGS)
l1qc_dct_c_LDFLAGS = $(LDFLAGS)
