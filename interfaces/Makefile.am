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

# $(info "MATHLIBS = $(MATHLIBS)")
# $(info "LDFLAGS = $(LDFLAGS)")
# $(info "MEXCPP = $(MEXCPPFLAGS)")
# $(info "MEXLIBS = $(MEXLDLIBS)")
# ----------------------------------------------------------------------#
# -------------------------------  ----------------------------------#

bin_PROGRAMS +=l1qc_dct_c
l1qc_dct_c_SOURCES = l1qc_dct_c.c libl1qc_dct.c
l1qc_dct_c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

l1qc_dct_c_CPPFLAGS =-I$(top_srcdir)/include  \
						 $(L1CCPPFLAGS)
l1qc_dct_c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
l1qc_dct_c_LDFLAGS = $(LDFLAGS)


lib_LTLIBRARIES += libl1qc_dct.la
libl1qc_dct_la_SOURCES	= libl1qc_dct.c
libl1qc_dct_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
							 $(L1CCPPFLAGS)
libl1qc_dct_la_CFLAGS	=  $(L1CCFLAGS) -fPIC \
							-std=c11 $(WARNFLAGS)
libl1qc_dct_la_LIBADD	= ../src/libl1c.la \
						  ../src/vcl/libvcl_math.la


libl1qc_dct_la_LDFLAGS	= -avoid-version -no-undefined $(LDFLAGS)\
							-Wl,--no-as-needed


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
l1qc_dct_mex_la_LDFLAGS  = -module -avoid-version -no-undefined -shrext $(MEXEXT) $(LDFLAGS)

l1qc_dct_mex_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
								$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_mex_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
								-std=c11 $(WARNFLAGS)

l1qc_dct_mex_la_LIBADD   =	libl1qc_dct.la \
							$(MEXLDLIBS)
endif


if ENABLE_PYTHON
# LIBEXT
# NOTES:
# https://www.gnu.org/software/libtool/manual/html_node/Link-mode.html
#
# To pass options to libtool, we use LDFLAGS.
#    -module makes it so we dont have to have the 'lib' prefix
#    -shrext .mexa64 lets us get a library without the .so, and with
#                    matlabs required extension.


lib_LTLIBRARIES += l1cPy.la

l1cPy_la_SOURCES  = l1cPy.c
l1cPy_la_LDFLAGS  = -module -avoid-version -no-undefined $(LDFLAGS) $(PYTHON_LIBS)

l1cPy_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
					 $(L1CCPPFLAGS) $(PYTHON_CPPFLAGS) $(NPY_CPPFLAGS)

l1cPy_la_CFLAGS   =  $(L1CCFLAGS) -fPIC \
						-std=c11 $(WARNFLAGS)

l1cPy_la_LIBADD   =	libl1qc_dct.la \
					../src/libl1c.la

endif


# Rules to copy the libs .libs into the current directory
# Idea from here:
# https://stackoverflow.com/questions/13948916/how-to-build-so-module-in-automake-and-a-libtool-using-project
# but modified because his regexp doesnt work.

.PHONY: convenience-link clean-convenience-link

# We need to ensure that we dont grep if lib_LTLTLIBRARIES is empty.
# Also, need to replace things like .so.0 with .so, so matlab can read it.
# lib_LTLLIBRARIES contains the .la files, which are just text descriptions.
convenience-link: $(lib_LTLIBRARIES)
	@echo "tlt libs $(lib_LTLIBRARIES)"
	@for soname in `test -z $$^ || $(EGREP) "^dlname=" $^|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
		new_soname=`echo $$soname|sed -e "s/so.*$$/so/g"`;                                         \
		if test ! -h "$$new_soname"; then                                                          \
		echo "$$soname: creating convenience link from $(abs_builddir)/.libs to $(abs_builddir)";  \
		echo "soname: $$soname";                                                                   \
		echo "new-soname: $$new_soname";                                                           \
		rm -f $(abs_builddir)/$$soname ;                                                           \
		test -e $(abs_builddir)/.libs/$$soname &&                                                  \
		$(LN_S) $(abs_builddir)/.libs/$$soname $$new_soname || echo "FAILED";                      \
		fi;                                                                                         \
	done


clean-convenience-link:
	rm -f *.so* *.mexa64 *.o *.lo *.la *.dll *.mexw64


all-local:: convenience-link
