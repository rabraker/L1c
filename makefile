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
MKLROOT        = /opt/intel/compilers_and_libraries_2019.1.144/linux/mkl
MATLAB         = matlab -nodesktop -nosplash -r
TEST_DATA_ROOT = $(CURDIR)/test_data
ML_BUILD_FCN   = L1qcTestData.build_all_test_data

ML_CMD         = "try;                                   \
                    ${ML_BUILD_FCN}('${TEST_DATA_ROOT}');\
                  catch E;                               \
                    fprintf('%s', E.message);            \
                    exit(1);                             \
                  end;                                   \
                 exit(0);"

SRC_DIR        = ./src
VCL_SRC_DIR    = ./src/vcl
TEST_SRC_DIR   = ./test
ML_INTERFACE   = ./interfaces

MEX_BUILD_DIR      = ./build/mex
MEX_LIB_DIR        = ./lib/mex
NATIVE_BUILD_DIR      = ./build
NATIVE_LIB_DIR        = ./lib

ifeq (${MAKECMDGOALS},mex)
BUILD_DIR      = $(MEX_BUILD_DIR)
LIB_DIR        = $(MEX_LIB_DIR)
else
BUILD_DIR      = $(NATIVE_BUILD_DIR)
LIB_DIR        = $(NATIVE_LIB_DIR)
endif
MEX_NAME       = $(ML_INTERFACE)/l1qc.mexa64
TEST_APP       = test_l1qc
ML_MEX         = $(ML_INTERFACE)/l1qc.mexa64

LIB_NAME       = $(LIB_DIR)/libl1qc.so
APP_ARCHIVE    = $(LIB_DIR)/libl1qc.a
VCL_ARCHIVE    = $(LIB_DIR)/libvcl_math.a

ifeq (${MAKECMDGOALS},test)
TEST_ARCHIVE   = test/lib/libtest_l1qc.a
else
TEST_ARCHIVE   =
endif


APP_SRC        = $(wildcard $(SRC_DIR)/*.c)
APP_CPP        = $(wildcard $(VCL_SRC_DIR)/*.cpp)
TEST_SRC       = $(wildcard $(TEST_SRC_DIR)/*.c)

APP_OBJ        = $(addprefix $(BUILD_DIR)/, $(notdir $(APP_SRC:.c=.o)))
VCL_OBJ        = $(addprefix $(BUILD_DIR)/, $(notdir $(APP_CPP:.cpp=.o)))
TEST_OBJ       = $(addprefix $(BUILD_DIR)/, $(notdir $(TEST_SRC:.c=.o)))
OBJ            = $(APP_OBJ) $(TEST_OBJ)

ifeq ($(DEBUG),1)
DBG          = -DDBUG -ggdb
OPT          =
VCL_OPT      =
else
DBG          =
OPT          = -O3 -msse3 -mavx2
VCL_OPT      = -O3 -mavx2 -mfma
endif
# check header files are in /usr/local/include
INCLUDE        = -I/usr/include                \
				 -I/usr/local/include          \
	           	 -Iinclude

# ---------- If making a mex file, add the needed includes and libs
ifeq ($(MAKECMDGOALS),mex)
MATLAB_CFLAGS  = -Wall -fPIC -D__MATLAB__ -fexceptions
MATLAB_INCLUDE = -I/usr/local/MATLAB/R2018b/extern/include
MATLAB_LDIR    = -L/usr/local/MATLAB/R2018b/bin/glnxa64 \

MATLAB_LIBS    = -lmex -lmat -lmx -lmwservices -lmwbuiltinsutil
endif

CFLAGS         =  $(OPT) $(DBG) -fPIC -Wall $(MAT_DEF)
CPP_VCL_FLAGS  =  -Iinclude/vcl $(DBG) $(VCL_OPT) -fabi-version=0 -fPIC

ifeq (${USE_MKL},1)
MATH_INCLUDE   =  -I${MKLROOT}/include
MATH_CFLAGS    =   -D_USEMKL_ -DMKL_ILP64 -m64
MATH_LDIR      =
MATH_ARCHIVE   =  ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
                  ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a  \
				  ${MKLROOT}/lib/intel64/libmkl_core.a
MATH_LIBS      = -lgomp -lpthread -lm -ldl

else
MATH_INCLUDE   =  -I/usr/local/OpenBlas/include
MATH_LDIR      =  -L/usr/local/OpenBlas/lib
MATH_LIBS      =  -lfftw3 -lopenblas -lm
endif

# put them all together
INCLUDE       +=  $(MATLAB_INCLUDE) $(MATH_INCLUDE)
CFLAGS        +=  $(MATLAB_CFLAGS) $(MATH_CFLAGS)

LIBDIR         = $(MATLAB_LDIR) $(MATH_LDIR)
LIBS           = $(MATH_LIBS)  $(MATLAB_LIBS) -lcheck -Wl,-rpath=/usr/local/lib
#N.B. ALL the archive files need to be enclosed in -Wl,--start-group -Wl,--end-group
# otherwise, linking may fail. See here:
# https://sourceware.org/binutils/docs/ld/Options.html
LDFLAGS       += -Wl,--start-group $(MATH_ARCHIVE) $(TEST_ARCHIVE) \
				$(APP_ARCHIVE) $(VCL_ARCHIVE) -Wl,--end-group \
				 $(LIBDIR) $(LIBS)


$(info $$DEF_MAT  [${ML_MEX}])


#-----------------------Targets ----------------------
.PHONY: clean help mex app_ar vcl_ar test_ar test_obj lib test test_data clean_all

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
	@echo "USE_MKL=1|0 (default 0). If 0, will use open_blas."
	@echo "DEBUG=1|0    (default 0). If 1, disables optimizations.\n"

mex: mex_obj
	${LD} -shared -o $(MEX_NAME) interfaces/l1qc.o ${LDFLAGS}

mex_obj:app_ar vcl_ar
	${CC} ${INCLUDE} ${CFLAGS} -shared -c -o interfaces/l1qc.o interfaces/l1qc.c

app_ar: $(APP_OBJ)
	$(AR) -rcs ${APP_ARCHIVE} $^

lib: $(APP_OBJ)
	@echo "building lib -------------"
	@echo "appobj = $(APP_OBJ)"
	$(LD) -o ${LIB_NAME} -shared $^ ${LIBDIR} $(LIBS)

test:test_ar app_ar vcl_ar
	$(LD) -o ${TEST_APP} $(LDFLAGS)

test_ar: $(TEST_OBJ)
	$(AR) -rcs $(TEST_ARCHIVE) $^

vcl_ar: $(VCL_OBJ)
	$(AR) -rcs $(VCL_ARCHIVE) $^

# $@ is to left side and $^ to the right of :
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c
	${CC} $(INCLUDE) ${CFLAGS} -c -o $@ $<

$(BUILD_DIR)/%.o: $(TEST_SRC_DIR)/%.c
	${CC} $(INCLUDE) ${CFLAGS} -c -o $@ $<

$(BUILD_DIR)/%.o: $(VCL_SRC_DIR)/%.cpp
	${CPP} $(CPP_VCL_FLAGS) -c -o $@ $<


test_data:
	${MATLAB} ${ML_CMD}

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



# ${CC_THREAD}
