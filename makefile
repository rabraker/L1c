# Unit-test Makefile

# Libraries have been installed in:
# /usr/local/lib
# If you ever happen to want to link against installed libraries
# in a given directory, LIBDIR, you must either use libtool, and
# specify the full pathname of the library, or use the `-LLIBDIR'
# flag during linking and do at least one of the following:
#    - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
#      during execution
#    - add LIBDIR to the `LD_RUN_PATH' environment variable
#      during linking
#    - use the `-Wl,-rpath -Wl,LIBDIR' linker flag
#    - have your system administrator add LIBDIR to `/etc/ld.so.conf'


# right now, need to use
# source /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl/bin/mklvars.sh intel64
# for mkl linking, the easiesy thing is to use:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/

# USE_MKL = false

CC             = gcc
CPP            = g++
LD             = gcc
AR             = ar

ifeq (${MAKECMDGOALS},mex)
DEF_MAT = -D__MATLAB__
else
DEF_MAT =
endif

# should check if this is set already
MEX_INSTALL_DIR = ../afm-cs/matlab-code/functions
MKLROOT        = /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl
# MKLROOT      = /opt/intel/compilers_and_libraries/linux/mkl
MATLAB_R       = /usr/local/MATLAB/R2018b
MATLAB         = matlab -nodesktop -nosplash -r
TEST_DATA_ROOT = $(CURDIR)/test_data
ML_BUILD_FCN   = L1qcTestData.build_all_test_data

ML_CMD         = "try;                                       \
                    ${ML_BUILD_FCN}('${TEST_DATA_ROOT}');    \
                  catch E;                                   \
                    fprintf('%s\n', E.message);              \
					for i=1:length(E.stack);                 \
						fprintf('Error in %s (line %d)\n',   \
						E.stack(i).name, E.stack(i).line);   \
					end;                                      \
                    exit(1);                                 \
                  end;                                       \
                 exit(0);"

SRC_DIR        = ./src
VCL_SRC_DIR    = ./src/vcl
TEST_SRC_DIR   = ./test
ML_INTERFACE   = ./interfaces

MEX_BUILD_DIR      = ./build/mex
MEX_LIB_DIR        = ./lib/mex
NATIVE_BUILD_DIR   = ./build
NATIVE_LIB_DIR     = ./lib

ifeq (${MAKECMDGOALS},mex)
BUILD_DIR      = $(MEX_BUILD_DIR)
APP_LIB_DIR    = $(MEX_LIB_DIR)
else
BUILD_DIR      = $(NATIVE_BUILD_DIR)
APP_LIB_DIR    = $(NATIVE_LIB_DIR)
endif
MEX_SRC        = $(ML_INTERFACE)/l1qc_dct.c
MEX_OBJ        = $(MEX_SRC:.c=.o)

# ML_MEX         = $(MEX_SRC:.c=.mexa64)
MEX_NAME       = $(MEX_SRC:.c=.mexa64)
TEST_APP       = test_l1qc

LIB_NAME       = $(APP_LIB_DIR)/libl1qc.so
APP_ARCHIVE    = $(APP_LIB_DIR)/libl1qc.a
VCL_ARCHIVE    = $(APP_LIB_DIR)/libvcl_math.a

ifeq (${MAKECMDGOALS},test)
TEST_LIB_DIR   = test/lib
TEST_ARCHIVE   = $(TEST_LIB_DIR)/libtest_l1qc.a
TEST_LIB       = -lcheck -lcjson
TEST_LDIR      = -Wl,-rpath=/usr/local/lib
else
TEST_ARCHIVE   =
TEST_LIB       =
TEST_LIB_DIR   =
TEST_LDIR   =
endif


APP_SRC        = $(wildcard $(SRC_DIR)/*.c)
APP_CPP        = $(wildcard $(VCL_SRC_DIR)/*.cpp)
TEST_SRC       = $(wildcard $(TEST_SRC_DIR)/*.c)

APP_OBJ        = $(addprefix $(BUILD_DIR)/, $(notdir $(APP_SRC:.c=.o)))
VCL_OBJ        = $(addprefix $(BUILD_DIR)/, $(notdir $(APP_CPP:.cpp=.o)))
TEST_OBJ       = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SRC:.c=.o)))
OBJ            = $(APP_OBJ) $(TEST_OBJ)

ifeq ($(DEBUG),1)
DBG            = -DDBUG -ggdb
OPT            =
VCL_OPT        =
else
DBG            =
# OPT          = -O3 -msse3 -mavx2
# VCL_OPT      = -O3 -mavx2 -mfma

# N.B. If we compile with unsupported instruction set flags, we will get
# a SIGILL error at runtime. For example, compiling on my desktop, which does not
# support avx2, only avx, with -maxv2 will cause a SIGILL. Thus, instead of specifying
# these flags explicitely, we can use -march=native -mtune=native. The result can be
# checked with
# gcc -dM -E -march=native -mtune=native - < /dev/null | egrep "AVX|SSE"
# which will print out a list of defines for enabled instructions.
OPT            = -O3 -march=native -mtune=native
VCL_OPT        = -O3 -march=native -mtune=native
endif
# check header files are in /usr/local/include
INCLUDE        = -I/usr/include                \
				 -I/usr/local/include          \
	           	 -Iinclude

# ---------- If making a mex file, add the needed includes and libs
ifeq ($(MAKECMDGOALS),mex)
MATLAB_CFLAGS  = -fPIC -D__MATLAB__ -fexceptions
MATLAB_INCLUDE = -I$(MATLAB_R)/extern/include
MATLAB_LDIR    = -L$(MATLAB_R)/bin/glnxa64 \

MATLAB_LIBS    = -lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
endif


# More warnings to look at:
# https://fastcompression.blogspot.com/2019/01/compiler-warnings.html
#  -Wstrict-aliasing -Wpointer-arith  \
#     -Wredundant-decls -Wmissing-prototypes\
#    -Wdeclaration-after-statement -Wc++-compat
#
# FAILS -Wredundant-decls -Wstrict-prototypes (seems like fault of mkl.h )
WARN_FLAGS     = -pedantic -Wall -Wextra -Wunused -Werror                   \
				 -Wcast-qual -Wcast-align -Winit-self -Wfloat-equal -Wundef \
				 -Wshadow -Wswitch-enum -Wvla

CFLAGS         =  $(OPT) $(DBG) -fopenmp -fPIC $(WARN_FLAGS) \
                  -std=c11 $(MAT_DEF)
CPP_VCL_FLAGS  =  -Iinclude/vcl $(DBG) $(VCL_OPT) $(WARN_FLAGS) -fabi-version=0 -fPIC

ifeq (${USE_OPENBLAS},1)
MATH_CFLAGS    = -D_HAVE_MM_MALLOC_
# MATH_CFLAGS    = -D_POSIX_C_SOURCE=200112L -D_HAVE_POSIX_MEMALIGN_
MATH_INCLUDE   =  -I/usr/local/OpenBlas/include
MATH_LDIR      =  -L/usr/local/OpenBlas/lib
MATH_LIBS      =  -lfftw3 -lgomp -lopenblas -lm

else
MATH_INCLUDE   =  -I${MKLROOT}/include
MATH_CFLAGS    =  -DMKL_ILP64  -D_USEMKL_  -m64
MATH_LDIR      =
# MATH_ARCHIVE   =  ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
# ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a  \
# ${MKLROOT}/lib/intel64/libmkl_core.a
# MATH_LIBS      = -lgomp -lpthread -lm -ldl
MATH_LIBS      = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_rt -lpthread -lm -ldl \
				-lfftw3 -lgomp
endif

# put them all together
INCLUDE       +=  $(MATLAB_INCLUDE) $(MATH_INCLUDE)
CFLAGS        +=   $(MATLAB_CFLAGS) $(MATH_CFLAGS)

LIBDIR         = $(MATLAB_LDIR) $(MATH_LDIR)  $(TEST_LDIR)
LIBS           = $(MATH_LIBS)  $(MATLAB_LIBS) $(TEST_LIB)
#N.B. ALL the archive files need to be enclosed in -Wl,--start-group -Wl,--end-group
# otherwise, linking may fail. See here:
# https://sourceware.org/binutils/docs/ld/Options.html
LDFLAGS       += -Wl,--start-group $(MATH_ARCHIVE) $(TEST_ARCHIVE) \
				$(APP_ARCHIVE) $(VCL_ARCHIVE) -Wl,--end-group \
				 $(LIBDIR) $(LIBS)


#-----------------------Targets ----------------------
.PHONY: clean help mex app_ar vcl_ar test_ar test_obj \
		test test_data clean_all install_mex

# test: $(TEST_APP)
help:
	@echo "Usage:"
	@echo ">>> make test [ARGS]"
	@echo "    Builds the test suite, which can be run by executing ./$(TEST_APP)"
	@echo ">>> make mex [ARGS]"
	@echo "    Builds the mex interface, which can be find in $(MEX_NAME)."
	@echo ">>> make test_data"
	@echo "    Builds test data required by the test suite. Running this command requires MATLAB."
	@echo ">>> make clean"
	@echo "    Removes all build artifacts, but leaves json test data and libraries."
	@echo ">>> make clean_all"
	@echo "    Removes all build artifacts including libraries, but leaves json test data."
	@echo ">>> make clean_pristine"
	@echo "    Removes all build artifacts including libraries and json test data."
	@echo "ARGS\n-----"
	@echo "USE_OPENBLAS=1|0 (default 0). If 1, will use open_blas."
	@echo "DEBUG=1|0    (default 0). If 1, disables optimizations.\n"

mex: $(MEX_NAME)

$(MEX_NAME):$(MEX_OBJ)
	${LD} -shared -o $(MEX_NAME) $< ${LDFLAGS}

$(MEX_OBJ):$(MEX_SRC) $(APP_ARCHIVE) $(VCL_ARCHIVE)
	${CC} ${INCLUDE} ${CFLAGS} -shared -c -o $@ $<

app_ar:$(APP_ARCHIVE)

vcl_ar:$(VCL_ARCHIVE)

app_lib: $(APP_OBJ)
	@echo "building lib -------------"
	$(LD) -o ${LIB_NAME} -shared $(LDFLAGS)

test:$(TEST_APP)


$(TEST_APP):$(TEST_ARCHIVE) $(APP_ARCHIVE) $(VCL_ARCHIVE)
	$(LD) -o ${TEST_APP} $(LDFLAGS)

$(APP_ARCHIVE): $(APP_OBJ) |$(APP_LIB_DIR)
	$(AR) -rcs ${APP_ARCHIVE} $^

$(TEST_ARCHIVE): $(TEST_OBJ) |$(TEST_LIB_DIR)
	$(AR) -rcs $(TEST_ARCHIVE) $^

$(VCL_ARCHIVE): $(VCL_OBJ) | $(APP_LIB_DIR)
	$(AR) -rcs $(VCL_ARCHIVE) $^

# $@ is to left side and $^ to the right of :
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c |$(BUILD_DIR)
	${CC} $(INCLUDE) ${CFLAGS} -c -o $@ $<

$(BUILD_DIR)/%.o: $(TEST_SRC_DIR)/%.c
	${CC} $(INCLUDE) ${CFLAGS} -c -o $@ $<

$(BUILD_DIR)/%.o: $(VCL_SRC_DIR)/%.cpp
	${CPP} $(CPP_VCL_FLAGS) -c -o $@ $<


$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(APP_LIB_DIR):
	@echo "app-lib-dir, test-lib-dir:$(TEST_LIB_DIR)"
	mkdir -p $(APP_LIB_DIR)

$(TEST_LIB_DIR):
	mkdir $(TEST_LIB_DIR)

test_data: |$(TEST_DATA_ROOT)
	${MATLAB} ${ML_CMD}

$(TEST_DATA_ROOT):
	mkdir $(TEST_DATA_ROOT)

install_mex:
# install the mex file
	cp $(MEX_NAME) $(MEX_INSTALL_DIR)/$(notdir $(MEX_NAME))
# install the mex-help file
	cp $(MEX_NAME:.mexa64=.m) $(MEX_INSTALL_DIR)/$(notdir $(MEX_NAME:.mexa64=.m))
# install the options builder
	cp interfaces/l1qc_opts.m $(MEX_INSTALL_DIR)/l1qc_dct_opts.m

clean:
	rm -f $(TEST_APP)
	rm -f $(MEX_BUILD_DIR)/*.o
	rm -f $(NATIVE_BUILD_DIR)/*.o
	rm -f $(ML_INTERFACE)/*.o
	rm -f $(ML_INTERFACE)/*.mexa64

clean_all:clean
	rm -f $(MEX_LIB_DIR)/*.*
	rm -f $(NATIVE_LIB_DIR)/*.*
	rm -f test/lib/*

clean_pristine:clean_all
	rm -f ${TEST_DATA_ROOT}/*.json


