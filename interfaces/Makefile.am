#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies

# MAT_CFLAGS = $(MKL_CFLAGS) $(OPT) $(DBG) -fopenmp -fPIC -Werror\
# 			-std=c11 -pedantic $(MAT_DEF)

MATLAB_R       = /usr/local/MATLAB/R2018b
# AM_CPPFLAGS = $(MKL_INCLUDE)

EXEEXT=.mexa64

bin_PROGRAMS=
$(info "MATHLIBS = $(MATHLIBS)")
$(info "LDFLAGS = $(LDFLAGS)")
$(info "MEXCPP = $(MEXCPPFLAGS)")


lib_LTLIBRARIES = libl1qc_dct.la
libl1qc_dct_la_SOURCES = libl1qc_dct.c
libl1qc_dct_la_CPPFLAGS = -I$(top_srcdir)/include  \
							$(L1CCPPFLAGS)
libl1qc_dct_la_CFLAGS =  $(L1CCFLAGS) -fPIC \
							-std=c11 $(WARNFLAGS)
libl1qc_dct_la_LIBADD = ../src/libl1c.la \
						../src/vcl/libvcl_math.la \
						$(MATHLIBS) -lpthread -lm -ldl -lgomp
libl1qc_dct_la_LDFLAGS = -Wl,--no-as-needed $(LDFLAGS)


if comp_with_mkl
if comp_with_mex
# LIBEXT
# NOTES:
# https://www.gnu.org/software/libtool/manual/html_node/Link-mode.html
#
# To pass options to libtool, we use LDFLAGS.
#    -module makes it so we dont have to have the 'lib' prefix
#    -shrext .mexa64 lets us get a library without the .so, and with
#                    matlabs required extension.



lib_LTLIBRARIES += l1qc_dct_mex.la

l1qc_dct_mex_la_SOURCES  = l1qc_dct_mex.c libl1qc_dct.c
l1qc_dct_mex_la_LDFLAGS = -module -avoid-version -shrext .mexa64 $(LDFLAGS)

l1qc_dct_mex_la_CPPFLAGS = -I$(top_srcdir)/include -I$(MATLAB_R)/extern/include \
								$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_mex_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
								-std=c11 $(WARNFLAGS)

l1qc_dct_mex_la_LIBADD    =	libl1qc_dct.la ../src/libl1c.la \
								../src/vcl/libvcl_math.la\
								$(MATHLIBS) -lpthread -lm -ldl -lgomp \
								-lmex -lmat -lmx -lmwservices -lmwbuiltinsutil

# bin_PROGRAMS += l1qc_dct_gnu
# l1qc_dct_gnu_SOURCES  = l1qc_dct_mex.c libl1qc_dct.c
# l1qc_dct_gnu_CPPFLAGS = -I$(top_srcdir)/include -I$(MATLAB_R)/extern/include \
# 							$(L1CCPPFLAGS) $(MEXCPPFLAGS)
# l1qc_dct_gnu_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
# 							-std=c11 $(WARNFLAGS)

# l1qc_dct_gnu_LDADD    =	../src/libl1c.a \
# 							../src/vcl/libvcl_math.a\
# 							$(MATHLIBS) -lpthread -lm -ldl -lgomp \
# 							-lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
# l1qc_dct_gnu_LDFLAGS = -shared $(LDFLAGS)

endif
endif

if comp_with_fftw
if comp_with_mex
bin_PROGRAMS += l1qc_dct_fftw

l1qc_dct_fftw_SOURCES = l1qc_fftw.c libl1qc_dct.c

l1qc_dct_fftw_CPPFLAGS = -I$(top_srcdir)/include -I$(MATLAB_R)/extern/include \
							$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_fftw_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS)  -fPIC \
						  -std=c11 $(WARNFLAGS)

l1qc_dct_fftw_LDADD =  	../src/libl1c.a \
						../src/vcl/libvcl_math.a\
						$(MATHLIBS) -lfftw3 -lpthread -lm -ldl -lgomp \
                         -lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
l1qc_dct_fftw_LDFLAGS = -shared -L$(MATLAB_R)/bin/glnxa64 \
						-Wl,--no-as-needed
endif
endif
