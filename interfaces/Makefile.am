#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies
clean-local:: clean-convenience-link
bin_PROGRAMS    =

# $(info "MATHLIBS = $(MATHLIBS)")
# $(info "LDFLAGS = $(LDFLAGS)")
# $(info "MEXCPP = $(MEXCPPFLAGS)")
# $(info "MEXLIBS = $(MEXLDLIBS)")
# -------------------------------------------------------------------#
# -------------------------------  ----------------------------------#

bin_PROGRAMS +=l1qc_dct_c
l1qc_dct_c_SOURCES = l1qc_dct_c.c
l1qc_dct_c_LDADD = ../src/libl1c.la ../src/vcl/libvcl_math.la

l1qc_dct_c_CPPFLAGS =-I$(top_srcdir)/include  \
						 $(L1CCPPFLAGS)
l1qc_dct_c_CFLAGS = $(L1CCFLAGS) -std=c11 $(WARNFLAGS)
l1qc_dct_c_LDFLAGS = $(LDFLAGS)


$(info "prefix= $(prefix)")
$(info "pyexecdir= $(pyexecdir)")
$(info "pyexecdir= $(PYEXECDIR)")

if comp_with_mex
# LIBEXT
# NOTES:
# https://www.gnu.org/software/libtool/manual/html_node/Link-mode.html
#
# To pass options to libtool, we use LDFLAGS.
#    -module makes it so we dont have to have the 'lib' prefix
#    -shrext .mexa64 lets us get a library without the .so, and with
#                    matlabs required extension.

mexdir = $(MEX_PREFIX)
mex_LTLIBRARIES = l1qc_dct.la \
				  breg_anistropic_TV.la

l1qc_dct_la_SOURCES  = l1qc_dct_mex.c
l1qc_dct_la_LDFLAGS  = -module -avoid-version -no-undefined -shrext $(MEXEXT) $(LDFLAGS)
l1qc_dct_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
								$(L1CCPPFLAGS) $(MEXCPPFLAGS)
l1qc_dct_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
								-std=c11 $(WARNFLAGS)
l1qc_dct_la_LIBADD   =	../src/libl1c.la \
							$(MEXLDLIBS)


breg_anistropic_TV_la_SOURCES  = breg_anistropic_TV_mex.c
breg_anistropic_TV_la_LDFLAGS  = -module -avoid-version -no-undefined -shrext $(MEXEXT) $(LDFLAGS)
breg_anistropic_TV_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
									$(L1CCPPFLAGS) $(MEXCPPFLAGS)
breg_anistropic_TV_la_CFLAGS   =  $(L1CCFLAGS) $(MEXCFLAGS) -fPIC \
									-std=c11 $(WARNFLAGS)
breg_anistropic_TV_la_LIBADD   = ../src/libl1c.la \
									$(MEXLDLIBS)


mex_DATA = l1qc_dct_opts.m l1qc_dct.m breg_anistropic_TV.m
endif


if ENABLE_PYTHON
# if $(pyexecdir) is empty, this will fail with
# libtool:   error: only absolute run-paths are allowed
python_PYTHON = __init__.py
pyexec_LTLIBRARIES = _l1cPy_module.la

_l1cPy_module_la_SOURCES  = l1cPy.c
_l1cPy_module_la_LDFLAGS  = -module -avoid-version $(LDFLAGS) $(PYTHON_LIBS)

_l1cPy_module_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
					 $(L1CCPPFLAGS) $(PYTHON_CPPFLAGS) $(NPY_CPPFLAGS)

_l1cPy_module_la_CFLAGS   =  $(L1CCFLAGS) -fPIC \
						-std=c11 $(WARNFLAGS)
_l1cPy_module_la_LIBADD   = $(top_builddir)/src/libl1c.la

endif



# Rules to copy the libs .libs into the current directory
# Idea from here:
# https://stackoverflow.com/questions/13948916/how-to-build-so-module-in-automake-and-a-libtool-using-project
# but modified because his regexp doesnt work.
all_libs = $(mex_LTLIBRARIES) $(pyexec_LTLIBRARIES)


.PHONY: convenience-link clean-convenience-link

# We need to ensure that we dont grep when lib_LTLTLIBRARIES is empty.
# This seems way overcomplicated...
# Also, need to replace things like .so.0 with .so, so matlab can read it.
# lib_LTLLIBRARIES contains the .la files, which are just text descriptions.
convenience-link: $(all_libs)
	@echo "----ltl_libs=$(all_libs)"
	if test `echo "x$(all_libs)"|sed s/\ //g` != x; then \
	for soname in `test -z $$^ || $(EGREP) "^dlname=" $^|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
		new_soname=`echo $$soname|sed -e "s/so.*$$/so/g"`;                                         \
		if test ! -h "$$new_soname"; then                                                          \
		echo "$$soname: creating convenience link from $(abs_builddir)/.libs to $(abs_builddir)";  \
		echo "soname: $$soname";                                                                   \
		echo "new-soname: $$new_soname";                                                           \
		rm -f $(abs_builddir)/$$soname ;                                                           \
		test -e $(abs_builddir)/.libs/$$soname &&                                                  \
		$(LN_S) $(abs_builddir)/.libs/$$soname $$new_soname || echo "FAILED";                      \
		fi;                                                                                         \
	done; \
	fi


clean-convenience-link:
	rm -f *.so* *.mexa64 *.o *.lo *.la *.dll *.mexw64 *.pyc


all-local:: convenience-link
