#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8
clean-local:: clean-convenience-link

# -------------------------------------------------------------------#
# -------------------------------  ----------------------------------#
example_data = $(abs_builddir)/examples/example_img_data.json

if ENABLE_PYTHON
# if $(pyexecdir) is empty, this will fail with
# libtool:   error: only absolute run-paths are allowed
python_PYTHON = __init__.py
pyexec_LTLIBRARIES = _l1cPy_module.la

_l1cPy_module_la_SOURCES  = l1cPy.c
_l1cPy_module_la_LDFLAGS  = -module -avoid-version $(LDFLAGS) $(PYTHON_LIBS)

_l1cPy_module_la_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include  \
							$(L1CCPPFLAGS) $(PYTHON_CPPFLAGS) \
							$(NPY_CPPFLAGS)
# -fPIC
_l1cPy_module_la_CFLAGS   =  $(L1CCFLAGS)
_l1cPy_module_la_LIBADD   = $(top_builddir)/src/libl1c.la

endif



.PHONY: clean-convenience-link copy-examples

# # We need to ensure that we dont grep when lib_LTLTLIBRARIES is empty.
# # This seems way overcomplicated...
# # Also, need to replace things like .so.0 with .so, so matlab can read it.
# # lib_LTLLIBRARIES contains the .la files, which are just text descrip



clean-convenience-link:
	rm -f *.so*  *.o *.lo *.la *.dll *.pyc examples/*


all-local::  copy-examples $(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT)) $(example_data)


$(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT)):.libs/$(pyexec_LTLIBRARIES)
	$(LN_S) -f $< $@


copy-examples:examples
	mkdir -p $(abs_builddir)/examples
	cp $(abs_srcdir)/examples/*.py $(abs_builddir)/examples/.
	$(LN_S) -f $(abs_builddir)/.libs/$(pyexec_LTLIBRARIES:.la=.$(SHRLIB_EXT))   $(abs_builddir)/examples/.


$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	echo "==************** $(example_data)"
	python $< $@
