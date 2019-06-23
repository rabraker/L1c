#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8



if ENABLE_PYTHON
# if $(pyexecdir) is empty, this will fail with
# libtool:   error: only absolute run-paths are allowed

# No wildcard in automake, because, portability.

mex_example_sources = cs20ng_grating.m   \
						mu_path_mask.m   \
						nesta_example.m  \
						play_l1qc_lb_and_TV_example.m

example_data = $(abs_top_builddir)/interfaces/mex/examples/example_img_data.json
example_build_dir = $(abs_top_builddir)/interfaces/mex/examples
mex_example_dest = $(addprefix $(example_build_dir)/, $(mex_example_sources))


EXTRA_DIST = $(mex_example_sources)

.PHONY: clean-convenience-link copy-examples


install-data-local:
	mkdir -p $(mexdir)/examples
	cp $(srcdir)/_examples/*.m $(mexdir)/examples/.



all-local:: $(example_build_dir) \
			$(example_data)      \
			$(py_example_dest)


# Create examples folder
$(example_build_dir):
	mkdir -p $@

# copy the mex examples scripts
$(example_build_dir)/%.m : $(srcdir)/%.m $(example_build_dir)
	$(LN_S) -f $< $@

# Create the example data.
$(example_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
	python3 $< $@


endif

clean-local:: clean-convenience-link

clean-convenience-link:
	rm -f $(example_build_dir)/*.m
