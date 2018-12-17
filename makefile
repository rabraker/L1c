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

APP_NAME       = test_l1magic

SRC_DIR        = src
BUILD_DIR      = build

SRC            = $(wildcard $(SRC_DIR)/*.c)
OBJ            = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC) )


DEBUG          = -DDBUG -ggdb
CC             = gcc
LD             = gcc
# check header files are in /usr/local/include
IFLAGS         = -I/usr/local/OpenBlas/include \
                 -I/usr/include                \
                 -I/usr/local/include          \
		          -Iinclude                     \


CFLAGS         = $(IFLAGS) $(DEBUG) -fPIC -msse3 -Wall

# For fftw: -lfftw3, #include <fftw3.h>
LDFLAGS        = -lfftw3 -lopenblas -lm -lcheck
LDFLAGS       += -Wl,-rpath=/usr/local/lib
LDFLAGS       +=  -L/usr/local/OpenBlas/lib


#-----------------------Targets ----------------------
.PHONY: clean all test_data clean_all

test:$(OBJ)
	$(LD) -o ${APP_NAME} $^  ${LDFLAGS}

all: test test_data
# $@ is to left side and $^ to the right of :
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c
	${CC} -g ${CFLAGS} -c -o $@ $<

test_data:
	${MATLAB} ${ML_CMD}

clean:
	rm -f $(APP_NAME)
	rm -f build/*.o

clean_all:clean
	rm -f ${TEST_DATA_ROOT}/*.json
