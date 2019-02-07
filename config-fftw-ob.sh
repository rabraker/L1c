#!/bin/bash

export CPPFLAGS="-I/usr/local/include -I/usr/local/MATLAB/R2018b/extern/include"
export LDFLAGS="-L/usr/local/lib -L/usr/local/MATLAB/R2018b/bin/glnxa64"

./configure --with-fftw3


# --with-mex
# ./configure $@

