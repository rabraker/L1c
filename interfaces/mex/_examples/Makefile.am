#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8

if ENABLE_PYTHON

# I want to be able to
# cd interfaces/mex/examples
# run them. I also want make clean to not destroy the example sources in
# an in-tree build. So we store the examples in _examples/, and symlink them to
# examples/

mex_example_sources = cs20ng_grating.m   \
						mu_path_mask.m   \
						nesta_example.m  \
						play_l1qc_lb_and_TV_example.m


example_data = $(abs_top_builddir)/interfaces/mex/examples/example_img_data.json
example_build_dir = $(abs_top_builddir)/interfaces/mex/examples
mex_example_dest = $(addprefix $(example_build_dir)/, $(mex_example_sources))


EXTRA_DIST = $(mex_example_sources)

.PHONY: clean-convenience-link copy-examples


install-data-local: $(MEX_PREFIX)/examples/l1c_mex_init_paths.m
	cp $(srcdir)/*.m $(MEX_PREFIX)/examples/.

# Use comma, rather than '/' for separator.
$(MEX_PREFIX)/examples/l1c_mex_init_paths.m:l1c_mex_init_paths.m.in $(MEX_PREFIX)/examples
	sed 's,\@MEX_ROOT\@,$(MEX_PREFIX),g' $< >$@

$(MEX_PREFIX)/examples:
	mkdir -p $(MEX_PREFIX)/examples


all-local::	$(example_data)      \
			$(mex_example_dest)


# copy the mex examples scripts. If we use ln -s, then l1c_mex_init_paths.m
# will always be newer than the (old) examples in $(srcdir).
$(example_build_dir)/%.m: $(srcdir)/%.m $(example_build_dir)/l1c_mex_init_paths.m
	cp -f $< $@

# Use comma, rather than '/' for separator.
$(example_build_dir)/l1c_mex_init_paths.m::$(srcdir)/l1c_mex_init_paths.m.in
	mkdir -p $(example_build_dir)
	sed 's,\@MEX_ROOT\@,$(abs_srcdir)/\.\.,g' $< >$@

# Create the example data.
$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	mkdir -p $(example_build_dir)
	python3 $< $@


clean-local::
	rm -rf $(example_build_dir)/*.m


endif
