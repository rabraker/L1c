# export LD_LIBRARY_PATH=/opt/python-versions/python3.4.2/lib/
#******************************************************************************
# gcc -o fgmpy.o fgm_python_interface.c -c -Wall -fPIC -I/usr/include -I/home/arnold/matlab/mex_sandbox  -I/usr/include/python3.4m

# gcc -shared -o fgmpy.so fgmpy.o

# Locating the root directory
# -g, debug
#CFLAGS = -Wall -I/usr/include -lgsl -lgslcblas -lm

ROOT=/home/arnold/matlab/afm-cs/c-src/l1magic
LIB_DIR = ${ROOT}
I_DIR = ${ROOT}

CC  = gcc -g
# CFLAGS = -c -Wall -fPIC -O -fexceptions
CFLAGS = -c
LFLAGS =  -L/usr/local/OpenBlas/lib  -lopenblas -lm
IFLAGS = -I/usr/local/OpenBlas/include -I${I_DIR}


COMMON_SRC = cgsolve.c
APP_NAME = cgsolve

# Rules for building the application and library



all: ${COMMON_SRC:.c=.o}
	$(CC)  ${COMMON_SRC:.c=.o} -o ${APP_NAME} ${LFLAGS}

# $@ is to left side and $^ to the right of :
${COMMON_SRC:.c=.o}:
	${CC} -o $@ ${@:.o=.c}  ${CFLAGS} ${IFLAGS}



debug:
	gdb $(APP_NAME)
clean:
	@rm -f *.o *.so *.mexa64



# ########################################################
# ROOT=/home/arnold/matlab/mex_sandbox
# # PY_ROOT = /home/arnold/.virtualenvs/rabraker-com-342/lib/python3.4

# CC  = gcc
# #CFLAGS = -Wall -I/usr/include -lgsl -lgslcblas -lm
# IFLAGS = -I/usr/include -I${ROOT}
# IFLAGS_MATLAB = ${IFLAGS} -I/usr/local/MATLAB/R2016b/extern/include

# CFLAGS = -c -Wall -fPIC
# EXTRAFLAGS_MATLAB = -O -D__MATLAB__ -fexceptions

# LFLAGS =  -lblas -lm

# LDIRS_MATLAB = -L/usr/local/MATLAB/R2016b/bin/glnxa64
# LFLAGS_MATLAB = ${LDIRS_MATLAB} ${LFLAGS} -lmex -lmat -lmx

# # APP_SRC  = fgm.c
# # APP_NAME = fgm
# APP_SRC = fgm_python_interface
# SRCS = fgm_python_interface utils
# APP_NAME = fgmpy
# MAT_EXT = mexa64

# # Rules for building the application and library
# #
# all: build
# matlab: build_matlab

# python: build_python

# build_matlab: toobj_matlab
# 	$(CC) -O2 -shared -o $(APP_NAME).${MAT_EXT} ${APP_SRC}.o ${LFLAGS_MATLAB}

# toobj_matlab:
# 	${CC} -O2 -o ${APP_NAME}.o ${APP_SRC} ${CFLAGS} ${IFLAGS_MATLAB} ${EXTRAFLAGS_MATLAB}

# run:
# 	./$(APP_NAME)
# debug:
# 	gdb $(APP_NAME)
# clean:
# 	@rm ${APP_NAME}.${EXT} ${APP_NAME}.o
