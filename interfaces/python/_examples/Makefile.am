#
#******************************************************************************
AUTOMAKE_OPTIONS = foreign 1.8


if ENABLE_PYTHON
# if $(pyexecdir) is empty, this will fail with
# libtool:   error: only absolute run-paths are allowed


example_data = $(abs_top_builddir)/interfaces/python/examples/example_img_data.json
example_build_dir = $(abs_top_builddir)/interfaces/python/examples
py_example_dest = $(addprefix $(example_build_dir)/, $(py_example_sources))

example_inst_dir = $(pkgpythondir)/examples
example_inst_data = $(example_inst_dir)/example_img_data.json

# No wildcard in automake, because, portability.
py_example_sources = dct_example.py           \
					run_TV_denoise_bregman.py \
					l1c_pyplot_utils.py

EXTRA_DIST = $(py_example_sources)

.PHONY: install-examples install-example-data


install-data-local:install-examples $(example_inst_data)

$(example_inst_data):$(example_inst_dir)
	python3 $(abs_top_srcdir)/test/build_CS20NG_example_data.py $@

install-examples:$(example_inst_dir)/l1c_py_init_path.py
	cp $(srcdir)/*.py $(pkgpythondir)/examples/.

$(example_inst_dir)/l1c_py_init_path.py:l1c_py_init_path.py.in $(example_inst_dir)
	$(SED) 's,\@L1C_PY_ROOT\@,$(pkgpythondir),g' $< >$@

$(example_inst_dir):
	$(MKDIR_P) $@



all-local:: $(example_data)      \
			$(py_example_dest)


# copy the python examples scripts
$(example_build_dir)/%.py : $(srcdir)/%.py $(example_build_dir) $(example_build_dir)/l1c_py_init_path.py
	cp -f $< $@

$(example_build_dir)/l1c_py_init_path.py:l1c_py_init_path.py.in $(example_build_dir)
	$(SED) 's,\@L1C_PY_ROOT\@,$(abs_builddir)/\.\./\.libs,g' $< >$@

$(example_build_dir):
	$(MKDIR_P) $@

# Create the example data.
$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py $(example_build_dir)
	python3 $< $@



clean-local::
	rm -f $(example_build_dir)/*.py \
			$(example_build_dir)/*.pyc \
			$(example_build_dir)/*.json

endif
