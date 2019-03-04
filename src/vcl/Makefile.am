
#******************************************************************************

# CPP_VCL_FLAGS  =  -Iinclude/vcl $(DBG) $(VCL_OPT) -fabi-version=0 -fPIC

# AM_CXXFLAGS = $(MKL_CFLAGS) $(OPT) $(DBG) -fopenmp -fPIC -Wall -Wextra -Wunused -Werror\
# 			-pedantic $(MAT_DEF)


# Using noinst_LTLLIBRARIES instead of noinst_LIBRARIES prevents the warning
# about linking against a static library when we compile and link the interfaces

noinst_LTLIBRARIES = libvcl_math.la
libvcl_math_la_SOURCES = vcl_math.cpp


# -O3 -march=native -mtune=native
# -fopenmp
AM_CPPFLAGS = -I$(top_srcdir)/include/vcl $(DBG) $(OPT) -fabi-version=0 -fPIC \
				-Wall -Wextra -Wunused -Werror\
	 			-pedantic
AM_CXXFLAGS = $(L1CCXXFLAGS)
