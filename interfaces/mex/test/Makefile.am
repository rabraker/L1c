#
#******************************************************************************

AUTOMAKE_OPTIONS = foreign 1.8
# clean-local:: clean-convenience-link
bin_PROGRAMS    =
AM_TESTS_ENVIRONMENT = 	ABS_TOP_SRCDIR='$(abs_top_srcdir)'; export ABS_TOP_SRCDIR;  \
						WITH_MEX='$(WITH_MEX)'; export WITH_MEX;

AM_LDFLAGS  = -module -avoid-version -no-undefined -shrext $(MEXEXT) $(LDFLAGS)
AM_CFLAGS   = $(L1CCFLAGS) $(MEXCFLAGS) -fPIC
AM_CPPFLAGS = -DBUILDING_DLL -I$(top_srcdir)/include \
	           -I$(top_srcdir)/interfaces/mex       \
				$(L1CCPPFLAGS) $(MEXCPPFLAGS)
LIBS        =  $(MEXLDLIBS) \
	$(top_builddir)/src/libl1c.la

# If we put l1c_mex_utils.la into LIBS above, we dont get dependency tracking. Ie,
# if we change l1c_mex_utils.c, it will get rebuilt, but our little test functions here
# wont get re-linked to it. Using LIBADD solves this, but is unfortunately more verbose.
mex_util_lalib = $(top_builddir)/interfaces/mex/l1c_mex_utils.la


if comp_with_mex
# LIBEXT
# NOTES:
# https://www.gnu.org/software/libtool/manual/html_node/Link-mode.html
#
# To pass options to libtool, we use LDFLAGS.
#    -module makes it so we dont have to have the 'lib' prefix
#    -shrext .mexa64 lets us get a library without the .so, and with
#                    matlabs required extension.


mexdir = $(MEX_PREFIX)/test
mex_LTLIBRARIES = TMU_mex_assert_num_inputs.la \
	TMU_mex_get_double_scalar_or_fail.la       \
	TMU_mex_get_double_array_or_fail.la        \
	TMU_mex_get_int_array_or_fail.la           \
	TMU_mex_assert_2Darray_with_size.la        \
	TMU_mex_assert_scalar_struct.la            \
	TMU_mex_get_double_from_struct_or_fail.la



TMU_mex_assert_num_inputs_la_SOURCES  = TMU_mex_assert_num_inputs.c ../l1c_mex_utils.h
TMU_mex_assert_num_inputs_la_LIBADD = $(mex_util_lalib)

TMU_mex_get_double_scalar_or_fail_la_SOURCES  = TMU_mex_get_double_scalar_or_fail.c ../l1c_mex_utils.h
TMU_mex_get_double_scalar_or_fail_la_LIBADD = $(mex_util_lalib)

TMU_mex_get_double_array_or_fail_la_SOURCES = TMU_mex_get_double_array_or_fail.c ../l1c_mex_utils.h
TMU_mex_get_double_array_or_fail_la_LIBADD =  $(mex_util_lalib)

TMU_mex_get_int_array_or_fail_la_SOURCES = TMU_mex_get_int_array_or_fail.c ../l1c_mex_utils.h
TMU_mex_get_int_array_or_fail_la_LIBADD = $(mex_util_lalib)

TMU_mex_assert_2Darray_with_size_la_SOURCES = TMU_mex_assert_2Darray_with_size.c ../l1c_mex_utils.h
TMU_mex_assert_2Darray_with_size_la_LIBADD = $(mex_util_lalib)

TMU_mex_assert_scalar_struct_la_SOURCES = TMU_mex_assert_scalar_struct.c ../l1c_mex_utils.h
TMU_mex_assert_scalar_struct_la_LIBADD = $(mex_util_lalib)

TMU_mex_get_double_from_struct_or_fail_la_SOURCES = TMU_mex_get_double_from_struct_or_fail.c ../l1c_mex_utils.h
TMU_mex_get_double_from_struct_or_fail_la_LIBADD = $(mex_util_lalib)


dist_check_SCRIPTS = l1c_assert_raises.m     \
             l1c_assert_error_free.m \
			 test_mex_interface.m    \
			 test_mex_utils.m        \
			 TMI_runner.m            \
			 clrs.m


endif
