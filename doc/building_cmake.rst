========
Building
========

.. _OpenBlas: https://github.com/xianyi/OpenBLAS
.. _FFTW3: http://fftw.org/
.. _check: https://github.com/libcheck/check
.. _cJSON: https://github.com/DaveGamble/cJSON

This project uses the CMake build system and Conan for some dependencies. 

To get started,

.. code-block:: bash

   git clone git@ggitlab.com:rabraker/L1c.git L1c && cd L1c

    sudo apt-get install libopenblas-openmp-dev python3-pip python3-dev \
        libfftw3-3 libfftw3-dev
    # create and activate python virtual environment if desired
    pip3 install -r tools/requirements.txt
    mkdir build && cd build
    conan install ..
    cmake ..
    make
    make test


Dependencies
==============
Running the full build requires the following dependencies

* OpenBlas_ This should be configured with openmp NOT pthread.
* FFTW3_
* check_ (used for unit tests). Tested with version 0.10 and 0.11
* cJSON_ (used for unit tests and example program). 
* GNU OpenMP (:code:`libgomp`) 
* Python3, Numpy, and Scipy are needed for unit tests and some examples.
* Python3-dev and python3-setuptools, if building the the python interface.
* Matplotlib is needed for plotting in some examples
* Valgrind (optional) to run the memory leak tests.




On FFTW and OpenBlas Threading
==============================
It is important that fftw and OpenBlas are build with OpenMP threading. If that is not the case, performance will greatly suffer. It's not entirely clear why, but
others have noticided similar issues, e.g.,

https://github.com/xianyi/OpenBLAS/issues/3187


mex bindings
==============
Mex bindings are not supported in the cmake build. See the deprecated Autools
build :doc:`building_autotools`.


Python bindings
===============

To build the python bindings, use :code:`--enable-python`:

.. code-block:: bash

   cmake -DBUILD_PYTHON_BINDINGS=ON ..

Building python bindings is supported for Python 3 (tested with 3.8). The proper compilation and linking flags as well as the installation location are obtained from the python3 on your path (via distutils.sysconfig). 


Unit Tests
==============
Almost all of the test data is generated in python and saved as json files in :code:`$(build_dir)/test/test_data/`.
To run the test suite, execute 

.. code-block:: bash

   make test

For more information about the tests, see :ref:`testing`.


