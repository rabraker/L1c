#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8
clean-local:: clean-convenience-link

# -------------------------------------------------------------------#
# -------------------------------  ----------------------------------#

if ENABLE_PYTHON
# if $(pyexecdir) is empty, this will fail with
# libtool:   error: only absolute run-paths are allowed
python_PYTHON = __init__.py
pyexec_LTLIBRARIES = _l1cPy_module.la

_l1cPy_module_la_SOURCES  = l1cPy.c
_l1cPy_module_la_LDFLAGS  = -module -avoid-version $(LDFLAGS) $(PYTHON_LIBS)

_l1cPy_module_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
					 $(L1CCPPFLAGS) $(PYTHON_CPPFLAGS) $(NPY_CPPFLAGS)

_l1cPy_module_la_CFLAGS   =  $(L1CCFLAGS) -fPIC
_l1cPy_module_la_LIBADD   = $(top_builddir)/src/libl1c.la

endif


# Rules to copy the libs .libs into the current directory
# Idea from here:
# https://stackoverflow.com/questions/13948916/how-to-build-so-module-in-automake-and-a-libtool-using-project
# but modified because his regexp doesnt work.
all_libs = $(pyexec_LTLIBRARIES)


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
	rm -f *.so*  *.o *.lo *.la *.dll *.pyc


all-local:: convenience-link
