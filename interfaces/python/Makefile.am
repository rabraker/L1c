#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8

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




.PHONY: clean-convenience-link copy-examples copy-lib

# # We need to ensure that we dont grep when lib_LTLTLIBRARIES is empty.
# # This seems way overcomplicated...
# # Also, need to replace things like .so.0 with .so, so matlab can read it.
# # lib_LTLLIBRARIES contains the .la files, which are just text descrip


# No wildcard in automake, because, portability.
py_example_sources = dct_example.py \
					run_TV_denoise_bregman.py

py_example_dest = $(addprefix $(abs_builddir)/examples/, $(py_example_sources))


all-local:: $(abs_builddir)/examples \
			copy-lib                 \
			$(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT)) \
			$(example_data)          \
			$(py_example_dest)


$(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT)):.libs/$(pyexec_LTLIBRARIES)
	$(LN_S) -f $< $@


$(abs_builddir)/examples:
	mkdir -p $@

$(abs_builddir)/examples/%.py:$(abs_srcdir)/examples/%.py $(abs_builddir)/examples
	$(LN_S) -f $< $@

copy-lib:
	$(LN_S) -f $(abs_builddir)/.libs/$(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT))   $(abs_builddir)/examples/.
	$(LN_S) -f $(abs_builddir)/.libs/$(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT))   $(abs_builddir)/.


$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	echo "==************** $(example_data)"
	python $< $@

endif

clean-local:: clean-convenience-link

clean-convenience-link:
	rm -f *.so*  *.o *.lo *.la *.dll *.pyc examples/*
