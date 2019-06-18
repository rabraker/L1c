#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8

# -------------------------------------------------------------------#
# -------------------------------  ----------------------------------#
example_data = $(abs_builddir)/examples/example_img_data.json
test_data = $(abs_builddir)/example_img_data.json

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

py_lib = $(pyexec_LTLIBRARIES:.la=$(SHRLIB_EXT))


all-local:: $(abs_builddir)/examples \
			copy-lib                 \
			$(example_data)          \
			$(test_data)          \
			$(py_example_dest)       \
			$(abs_builddir)/test_L1cPy.py \
			$(abs_builddir)/examples/$(py_lib) \
			$(abs_builddir)/$(py_lib)

# Ensure the .la library exists .libs
.libs/$(pyexec_LTLIBRARIES):$(pyexec_LTLIBRARIES)


# copy lib to this folder, for tests.
$(abs_builddir)/$(py_lib):.libs/$(pyexec_LTLIBRARIES)
	$(LN_S) -f $(abs_builddir)/$(<:.la=$(SHRLIB_EXT)) $@


# copy lib to examples folder
$(abs_builddir)/examples/$(py_lib):.libs/$(pyexec_LTLIBRARIES)
	$(LN_S) -f $(abs_builddir)/$(<:.la=$(SHRLIB_EXT)) $@

# Copy the test script. Python path gets screwed up when using a symlink,
# so just copy it.
$(abs_builddir)/test_L1cPy.py:$(abs_srcdir)/test_L1cPy.py
	cp -f $< $@

# Create examples folder
$(abs_builddir)/examples:
	mkdir -p $@

# copy the python examples scripts
$(abs_builddir)/examples/%.py:$(abs_srcdir)/_examples/%.py $(abs_builddir)/examples
	$(LN_S) -f $< $@

# Create the example data.
$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	python3 $< $@

$(test_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	python3 $< $@


endif

clean-local:: clean-convenience-link

clean-convenience-link:
	rm -rf *.so*  *.o *.lo *.la *.dll *.pyc examples/
