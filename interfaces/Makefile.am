#
#******************************************************************************
# If configure is run with --with-debug, then
#        OPT = -O0 -ggdb and OPT = ''.
# Otherwise,
#       DBG = '' and OPT = -O3 -march=native -mtune=native

AUTOMAKE_OPTIONS = foreign 1.8 no-dependencies

MAT_CFLAGS = $(MKL_CFLAGS) $(OPT) $(DBG) -fopenmp -fPIC -Wall -Wextra -Wunused -Werror\
			-std=c11 -pedantic $(MAT_DEF)

# MKLROOT        = /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl
MATLAB_R       = /usr/local/MATLAB/R2018b
AM_CPPFLAGS = $(MKL_INCLUDE) -I$(top_srcdir)/include -I$(MATLAB_R)/extern/include $(L1C_FLAGS)

$(info $(L1C_FLAGS))
AM_CFLAGS = $L1C_FLAGS
bin_PROGRAMS = l1qc_dct_gnu
l1qc_dct_gnu_SOURCES = l1qc_dct.c
l1qc_dct_gnu_CFLAGS  = $(L1C_CFLAGS) $(OPT) $(DBG) $(CFLAGS) -fopenmp -fPIC \
						-std=c11 -pedantic $(MAT_DEF)

l1qc_dct_gnu_LDADD =  		../src/libl1c.a \
							../src/vcl/libvcl_math.a\
						-lmkl_rt -lpthread -lm -ldl -lgomp -lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
l1qc_dct_gnu_LDFLAGS = -shared -L$(MATLAB_R)/bin/glnxa64 \
						 -Wl,--no-as-needed
