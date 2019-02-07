#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
clean-local:: clean-convenience-link

bin_PROGRAMS    =
lib_LTLIBRARIES =

$(info "MATHLIBS = $(MATHLIBS)")
$(info "LDFLAGS = $(LDFLAGS)")
$(info "MEXCPP = $(MEXCPPFLAGS)")

# ----------------------------------------------------------------------#
# ------------------------------- MKL ----------------------------------#
if comp_with_mkl

bin_PROGRAMS +=l1qc_dct_c
l1qc_dct_c_SOURCES = l1qc_dct_c.c libl1qc_dct.c
l1qc_dct_c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la \
					$(MATHLIBS)  -lm -ldl -lgomp
# -lpthread
l1qc_dct_c_CPPFLAGS =-I$(top_srcdir)/include  \
					-D_AX_MKL_ $(L1CCPPFLAGS)
l1qc_dct_c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
l1qc_dct_c_LDFLAGS = $(LDFLAGS) -pthread



lib_LTLIBRARIES += libl1qc_dct_mkl.la
libl1qc_dct_mkl_la_SOURCES  = libl1qc_dct.c
libl1qc_dct_mkl_la_CPPFLAGS = -I$(top_srcdir)/include  \
							-D_AX_MKL_ $(L1CCPPFLAGS)
libl1qc_dct_mkl_la_CFLAGS   =  $(L1CCFLAGS) -fPIC \
							-std=c11 $(WARNFLAGS)
libl1qc_dct_mkl_la_LIBADD   = ../src/libl1c.la \
							../src/vcl/libvcl_math.la \
							$(MATHLIBS)  -lm -ldl -lgomp
libl1qc_dct_mkl_la_LDFLAGS  = -module 	 $(LDFLAGS) -pthread

# -lpthread

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

l1qc_dct_mex_la_SOURCES  = l1qc_dct_mex.c
l1qc_dct_mex_la_LDFLAGS  = -module -avoid-version -shrext .mexa64 $(LDFLAGS) -pthread

l1qc_dct_mex_la_CPPFLAGS = -I$(top_srcdir)/include  \
								$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_mex_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
								-std=c11 $(WARNFLAGS)

l1qc_dct_mex_la_LIBADD   =	libl1qc_dct_mkl.la \
							$(MATHLIBS)  -lm -ldl -lgomp \
							-lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
# -lpthread
endif
endif

# ----------------------------------------------------------------------#
# --------------------------- FFTW/OpenBlas ----------------------------#
if comp_with_fftw

bin_PROGRAMS +=l1qc_dct_c_fftw
l1qc_dct_c_fftw_SOURCES = l1qc_dct_c.c libl1qc_dct.c
l1qc_dct_c_fftw_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la \
							$(MATHLIBS)  -lm -ldl -lgomp
# -lpthread
l1qc_dct_c_fftw_CPPFLAGS =-I$(top_srcdir)/include  \
							-D_AX_FFTW_ $(L1CCPPFLAGS)
l1qc_dct_c_fftw_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
l1qc_dct_c_fftw_LDFLAGS = $(LDFLAGS) -pthread


lib_LTLIBRARIES += libl1qc_dct_fftw.la
libl1qc_dct_fftw_la_SOURCES  = libl1qc_dct.c
libl1qc_dct_fftw_la_CPPFLAGS = -I$(top_srcdir)/include  \
								-D_AX_FFTW_ $(L1CCPPFLAGS)
libl1qc_dct_fftw_la_CFLAGS   =  $(L1CCFLAGS) -fPIC \
									-std=c11 $(WARNFLAGS)
libl1qc_dct_fftw_la_LIBADD   = ../src/libl1c.la \
								../src/vcl/libvcl_math.la \
								$(MATHLIBS)  -lm -ldl -lgomp
libl1qc_dct_fftw_la_LDFLAGS  = -module 	 $(LDFLAGS) -pthread


if comp_with_mex
bin_PROGRAMS += l1qc_dct_fftw

l1qc_dct_fftw_SOURCES  = l1qc_fftw.c libl1qc_dct.c

l1qc_dct_fftw_CPPFLAGS = -I$(top_srcdir)/include  \
							-D_AX_FFTW_ $(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_fftw_CFLAGS   = $(L1CCFLAGS) $(MEXCFLAGS)  -fPIC \
						  -std=c11 $(WARNFLAGS)

l1qc_dct_fftw_LDADD    = ../src/libl1c.a \
							../src/vcl/libvcl_math.a\
							$(MATHLIBS) -lfftw3 -lpthread -lm -ldl -lgomp \
							-lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
l1qc_dct_fftw_LDFLAGS  = -shared -L$(MATLAB_R)/bin/glnxa64 \
							-Wl,--no-as-needed
endif
endif


# Rules to copy the libs .libs into the current directory
# Idea from here:
# https://stackoverflow.com/questions/13948916/how-to-build-so-module-in-automake-and-a-libtool-using-project
# but modified because his regexp doesnt work.

.PHONY: convenience-link clean-convenience-link

convenience-link: $(lib_LTLIBRARIES)
	@echo "tlt libs $(lib_LTLIBRARIES)"
	@for soname in `$(EGREP) "^dlname=" $^|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
		echo "$$soname: creating convenience link from $(abs_builddir)/.libs to $(abs_builddir)"; \
		echo "soname: $$soname";\
		rm -f $(abs_builddir)/$$soname ; \
		test -e $(abs_builddir)/.libs/$$soname && \
		$(LN_S) $(abs_builddir)/.libs/$$soname $$soname || echo "FAILED";\
	done

# Im not sure what exactly this was supposed to do, but it doesnt work.
# @for soname in `echo | $(EGREP) "^dlname=" $^ | $(SED) -e "s|^dlname='\(.*\)'|\1|"`; do  \
#   	echo "$$soname: creating convenience link from $(abs_builddir)/.libs to $(top_builddir)"; \
# 	rm -f $(top_builddir)/$$soname ; \
# 	test -e $(abs_builddir)/.libs/$$soname && \
# 	cd $(top_builddir) && \
# 	$(LN_S) $(abs_builddir)/.libs/$$soname $$soname || true;\
# done

clean-convenience-link:
	rm -f *.so* *.mexa64

# @echo "input to clean-conf:  $(lib_LTLIBRARIES)"
# @for soname in `$(EGREP) "^dlname=" $(lib_LTLIBRARIES)|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
# 	echo "$$soname: cleaning convenience links"; \
# done

# 		test -L $(abs_builddir)/$$soname && rm -f $(abs_builddir)/$$soname || true; \
#
all-local:: convenience-link
