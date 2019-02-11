#!/bin/bash

# if test x$HOSTNAME = xzeus;
# then
    export CPPFLAGS="-I/usr/local/include -I/usr/local/MATLAB/R2018b/extern/include"
    export LDFLAGS="-L/usr/local/lib -L/usr/local/MATLAB/R2018b/bin/glnxa64"
# elif test x$HOSTNAME = xegoipse;
# then
#     export CPPFLAGS="-I/usr/local/OpenBlas/include -I/usr/local/MATLAB/R2018b/extern/include"
#     export LDFLAGS="-L/usr/local/OpenBlas/lib -L/usr/local/MATLAB/R2018b/bin/glnxa64"
# fi

./configure --with-fftw3 --with-mex "${@}"


# --with-mex
# ./configure $@

