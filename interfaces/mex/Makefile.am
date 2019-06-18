#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8
clean-local:: clean-convenience-link

bin_PROGRAMS    =

AM_LDFLAGS      = -module -avoid-version -no-undefined -shrext .$(MEXEXT) $(LDFLAGS)
AM_CPPFLAGS     = -DBUILDING_DLL -I$(top_srcdir)/include -I./ -I../ \
					$(L1CCPPFLAGS) $(MEXCPPFLAGS)
#  -fPIC
AM_CFLAGS       =  $(L1CCFLAGS) $(MEXCFLAGS)
LIBS            =  $(top_builddir)/src/libl1c.la \
							$(MEXLDLIBS)
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
				  breg_anistropic_TV.la \
				  nesta_dctTV.la

noinst_LTLIBRARIES = l1c_mex_utils.la

l1c_mex_utils_la_SOURCES = l1c_mex_utils.c l1c_mex_utils.h

l1qc_dct_la_SOURCES  = l1qc_dct_mex.c
l1qc_dct_la_LIBADD = l1c_mex_utils.la

breg_anistropic_TV_la_SOURCES  = breg_anistropic_TV_mex.c
breg_anistropic_TV_la_LIBADD = l1c_mex_utils.la

nesta_dctTV_la_SOURCES  = nesta_dctTV_mex.c
nesta_dctTV_la_LIBADD = l1c_mex_utils.la

dist_mex_DATA = l1qc_dct_opts.m \
				l1qc_dct.m     \
				breg_anistropic_TV.m \
				nesta_opts.m



# Rules to copy the libs .libs into the current directory
# Idea from here:
# https://stackoverflow.com/questions/13948916/how-to-build-so-module-in-automake-and-a-libtool-using-project
# but modified because his regexp doesnt work.
all_libs = $(mex_LTLIBRARIES)

SUFFIXES=.$(MEXEXT)

.PHONY: convenience-link clean-convenience-link copy-examples

# We need to ensure that we dont grep when lib_LTLTLIBRARIES is empty.
# This seems way overcomplicated...
# Also, need to replace things like .so.0 with .so, so matlab can read it.
# lib_LTLLIBRARIES contains the .la files, which are just text descriptions.
# convenience-link: $(all_libs)
# 	@echo "----ltl_libs=$(all_libs)"
# 	if test `echo "x$(all_libs)"|sed s/\ //g` != x; then \
# 	for soname in `test -z $$^ || $(EGREP) "^dlname=" $^|sed -e "s/^.*=//g"|sed -e "s/'//g"`; do  \
# 		new_soname=`echo $$soname|sed -e "s/so.*$$/so/g"`;                                         \
# 		if test ! -h "$$new_soname"; then                                                          \
# 		echo "$$soname: creating convenience link from $(abs_builddir)/.libs to $(abs_builddir)";  \
# 		echo "soname: $$soname";                                                                   \
# 		echo "new-soname: $$new_soname";                                                           \
# 		rm -f $(abs_builddir)/$$soname ;                                                           \
# 		test -e $(abs_builddir)/.libs/$$soname &&                                                  \
# 		$(LN_S) $(abs_builddir)/.libs/$$soname $$new_soname || echo "FAILED";                      \
# 		fi;                                                                                         \
# 	done; \
# 	fi

mex_example_sources = cs20ng_grating.m  \
						mu_path_mask.m  \
						nesta_example.m  \
						play_l1qc_lb_and_TV_example.m

mex_example_dest := $(addprefix $(abs_builddir)/examples/, $(mex_example_sources))

mex_libs = $(mex_LTLIBRARIES:.la=.$(MEXEXT))

mex_lib_example_dest :=$(addprefix $(abs_builddir)/examples/, $(mex_libs))
mex_lib_dest := $(addprefix $(abs_builddir)/, $(mex_libs))

$(info "----------mex dest: $(mex_lib_dest)")


all-local:: $(abs_builddir)/examples \
			$(mex_example_dest)      \
			$(mex_lib_dest)           \
			$(mex_lib_example_dest)



$(abs_builddir)/%.$(MEXEXT):%.la
	$(LN_S) -f $(abs_builddir)/.libs/$(<:.la=.mexa64) $@

$(abs_builddir)/examples/%.$(MEXEXT):%.la
	$(LN_S) -f $(abs_builddir)/.libs/$(<:.la=.mexa64) $@


# Link (or copy) the example scripts to the examples folder in the build dir.
$(abs_builddir)/examples/%.m:$(abs_srcdir)/examples/%.m $(abs_builddir)/examples
	$(LN_S) -f $< $@

$(abs_builddir)/examples:
	mkdir -p $@

endif


clean-convenience-link:
	rm -f *.so* *.mexa64 *.o *.lo *.la *.dll *.mexw64 *.pyc examples/*.$(MEXEXT)
