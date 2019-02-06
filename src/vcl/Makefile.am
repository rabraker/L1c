
#******************************************************************************

# CPP_VCL_FLAGS  =  -Iinclude/vcl $(DBG) $(VCL_OPT) -fabi-version=0 -fPIC

# AM_CXXFLAGS = $(MKL_CFLAGS) $(OPT) $(DBG) -fopenmp -fPIC -Wall -Wextra -Wunused -Werror\
# 			-pedantic $(MAT_DEF)


# lib_LTLIBRARIES = libvcl_math.la
# libvcl_math_la_SOURCES = vcl_math.cpp

noinst_LTLIBRARIES = libvcl_math.la
libvcl_math_la_SOURCES = vcl_math.cpp


# -O3 -march=native -mtune=native
AM_CPPFLAGS = -I$(top_srcdir)/include/vcl $(DBG) $(OPT) -fopenmp -fabi-version=0 -fPIC \
				-Wall -Wextra -Wunused -Werror\
	 			-pedantic
AM_CXXFLAGS = $(L1CCXXFLAGS)
