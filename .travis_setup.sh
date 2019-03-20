#!/bin/bash

sudo apt-get install -y libopenblas-base libopenblas-dev \
	python3-numpy python3-scipy libfftw3-3 libfftw3-dev


check_release=0.12.0
wget "https://github.com/libcheck/check/releases/download/${check_release}/check-${check_release}.tar.gz"
tar -xvvf "check-${check_release}.tar.gz"
pushd $(pwd)
cd "check-${check_release}"
./configure
make
sudo make install
popd

git clone https://github.com/DaveGamble/cJSON.git
pushd $(pwd)
cd cJSON
make
sudo make install
popd

sudo ldconfig
