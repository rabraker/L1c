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
$(info "MEXLIBS = $(MEXLDLIBS)")
# ----------------------------------------------------------------------#
# -------------------------------  ----------------------------------#
if compile_dct_example
# -lm -ldl -lgomp
bin_PROGRAMS +=l1qc_dct_c
l1qc_dct_c_SOURCES = l1qc_dct_c.c libl1qc_dct.c
l1qc_dct_c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

l1qc_dct_c_CPPFLAGS =-I$(top_srcdir)/include  \
						 $(L1CCPPFLAGS)
l1qc_dct_c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
l1qc_dct_c_LDFLAGS = $(LDFLAGS)


lib_LTLIBRARIES += libl1qc_dct.la
libl1qc_dct_la_SOURCES	= libl1qc_dct.c
libl1qc_dct_la_CPPFLAGS = -I$(top_srcdir)/include  \
							 $(L1CCPPFLAGS)
libl1qc_dct_la_CFLAGS	=  $(L1CCFLAGS) -fPIC \
						-std=c11 $(WARNFLAGS)
libl1qc_dct_la_LIBADD	= ../src/libl1c.la \
						../src/vcl/libvcl_math.la


libl1qc_dct_la_LDFLAGS	= -module	 $(LDFLAGS)\
	-Wl,--no-as-needed -Wl,--rpath=/opt/intel/compilers_and_libraries_2019.1.144/linux/mkl/lib/intel64_lin

# This worked: note libmkl_avx.so
# gcc -shared  -fPIC -DPIC  .libs/libl1qc_dct_la-libl1qc_dct.o  -Wl,--whole-archive ../src/.libs/libl1c.a ../src/vcl/.libs/libvcl_math.a -Wl,--no-whole-archive  -L/opt/intel/mkl/lib/intel64 /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.so /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.so /opt/intel/mkl/lib/intel64/libmkl_core.so /opt/intel/mkl/lib/intel64/libmkl_avx.so -lgomp -lm -ldl  -lcjson -O3 -march=native -mtune=native -fopenmp -m64 -Wl,--no-as-needed -Wl,--rpath,/opt/intel/mkl/lib/intel64   -fopenmp -Wl,-soname -Wl,libl1qc_dct.so.0 -o libl1qc_dct.so.0

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
l1qc_dct_mex_la_LDFLAGS  = -module -avoid-version -shrext $(MEXEXT) $(LDFLAGS)

l1qc_dct_mex_la_CPPFLAGS = -I$(top_srcdir)/include  \
								$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_mex_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
								-std=c11 $(WARNFLAGS)

l1qc_dct_mex_la_LIBADD   =	libl1qc_dct.la \
							$(MEXLDLIBS)
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
	rm -f *.so* *.mexa64 *.o *.lo *.la

# @echo "input to clean-conf:  $(lib_LTLIBRARIES)"
# @for soname in `$(EGREP) "^dlname=" $(lib_LTLIBRARIES)|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
# 	echo "$$soname: cleaning convenience links"; \
# done

# 		test -L $(abs_builddir)/$$soname && rm -f $(abs_builddir)/$$soname || true; \
#
all-local:: convenience-link
