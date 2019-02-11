#!/bin/bash

export CPPFLAGS="-I/opt/intel/mkl/include"
export LDFLAGS="-L/opt/intel/mkl/lib/intel64"

./configure --with-mkl $@


