#!/bin/bash

export CPPFLAGS="-I/opt/intel/mkl/include -I/usr/local/R2018b/extern/include" 
export LDFLAGS="-L/opt/intel/mkl/lib/intel64 -L/usr/local/MATLAB/R2018b/bin/glnxa64" 

./configure $@


