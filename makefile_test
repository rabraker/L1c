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

APP_NAME      = test_l1magic

SRC_DIR   = src
BUILD_DIR = build

SRC       = $(wildcard $(SRC_DIR)/*.c)
OBJ       = $(patsubst $(SRC_DIR)/%.c, $(BUILD_DIR)/%.o, $(SRC) )


DEBUG     = -DDBUG -ggdb
CC        = gcc
LD        = gcc
# check header files are in /usr/local/include
IFLAGS    = -I/usr/local/OpenBlas/include \
            -I/usr/include                \
            -I/usr/local/include          \
		    -Iinclude                     \


CFLAGS    = $(IFLAGS) $(DEBUG) -fPIC -msse3 -Wall

# For fftw: -lfftw3, #include <fftw3.h>
LDFLAGS   = -lfftw3 -lopenblas -lm -lcheck
LDFLAGS   += -Wl,-rpath=/usr/local/lib
LDFLAGS   +=  -L/usr/local/OpenBlas/lib


#---------------------------------------------------------- Targets
.PHONY: clean all test

all:$(OBJ)
	$(LD) -o ${APP_NAME} $^  ${LDFLAGS}

# $@ is to left side and $^ to the right of :
$(BUILD_DIR)/%.o : $(SRC_DIR)/%.c
	${CC} -g ${CFLAGS} -c -o $@ $<


clean:
	rm -f $(APP_NAME)
	rm -f build/*.o
