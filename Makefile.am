
AUTOMAKE_OPTIONS = foreign

# The order of these is important: e.g., we need to run src/vcl/Makefile
# before we run src/Makefile
SUBDIRS = src/vcl             \
			src               \
			interfaces/python/_examples \
		  	interfaces/python \
			interfaces/mex    \
			interfaces/mex/_examples    \
            interfaces/mex/test \
			examples/c        \
			test              \
			include           \
			include/vcl



# libtool requires some M4 scripts to be added to the source tree.  Make sure that
# Autoconf knows where to find them.
ACLOCAL_AMFLAGS = -I m4
