========
Building
========

.. _OpenBlas: https://github.com/xianyi/OpenBLAS
.. _ATLAS: http://math-atlas.sourceforge.net
.. _FFTW3: http://fftw.org/
.. _check: https://github.com/libcheck/check
.. _cJSON: https://github.com/DaveGamble/cJSON

This project uses the GNU Autotools build system. Considering that we are still at "version 0.0", there is not a release tarball yet, so you will need Autotools installed. I typically test on Debian Linux. I have also succesfully built the project on Windows 7 with MinGW64.

The four steps are 

.. code-block:: bash

   git clone git@ggitlab.com:rabraker/L1c.git L1c && cd L1c
   ./bootsrap                   # requires autotools
   ./configure [--with-fftw3 [--with-fftw3-threads=[ARG]]] [--with-mkl] [--with-mex]
   make && make install

Running :code:`./configure --help` will give you a list of options. Depending on where the dependencies were installed, it may be necessary to augment the :code:`./configure` step with the right environment variables. 


Dependencies
==============
Running the full build requires the following dependencies

* OpenBlas_ or ATLAS_
* FFTW3_
* check_ (used for unit tests). Tested with version 0.10 and 0.11
* cJSON_ (used for unit tests and example program). 
* GNU OpenMP (:code:`libgomp`) 
* Python3, Numpy, and Scipy are needed for unit tests and some examples.
* Python3-dev and python3-setuptools, if building the the python interface.
* Matplotlib is needed for plotting in some examples
* Valgrind (optional) to run the memory leak tests.

Details and suggestions follow. Or run :code:`./configure --help` if this is old hat to you.

Presently, the configure script will not finish if it does not find :code:`cJSON` or :code:`check`.

Core Code
==============
The core code with examples requires FFTW and either OpenBlas or ATLAS. 
In my experiments, if OpenBlas/Atlas and FFTW3 are compiled with support for threading and the avx and sse instruction sets, these combinations have slightly better performance than MKL. It is likely desirable (certainly for ATLAS) to compile these yourself, since the libraries availible through your distribution may be older or not have been compiled with full optimization. For ATLAS, the compilation optimizes the binary for your specific computer, so it makes no sense to obtain a library built by somebody else. 

By default, the configure script will search for first OpenBlas, then ATLAS. If you have both installed, you can instruct :code:`configure` to choose one over the other by setting the envirnmental variable :code:`BLAS_LIB=[openblas][satlas]`. The key requirement for the BLAS library is that it has the extension :code:`cblas_daxpby`, or in the case of ATLAS, :code:`catlas_daxpby`. Thus, in principle, you can use another BLAS library with this extension, though that is untested.

FFTW3 threading comes in three flavors: (1) OpenMP threading (`libfftw_omp`), (2) POSIX threads (`libfftw3_threads`) and (3), as a single, combined library. To choose a particular threading version, use the flag :code:`--with-fftw3-threads=[combined][omp][threads][yes][no]`. The :code:`yes` option will search over all three possibilities. On Windows with MinGW, FFTW3 only compiles threading into a combined library.

To enable compiling a Matlab mex functions, use :code:`--with-mex`. This requires that the Mathworks  script :code:`mexext` is on your path, or that you export the environmental variable :code:`MEXEXT` with the proper mex extension for your platform.

In general, it is likely (almost certain for mex) that these libraries will not be on your default search path, so you must tell :code:`configure` where they are by defining them in the :code:`CPPFLAGS` and :code:`LDFLAGS` environmental variables. See the example below. 

Example Builds
==============
On my machine, :code:`FFTW3`, :code:`cJSON`, and :code:`check` are installed in :code:`/usr/local`. ATLAS is installed in :code:`/usr/local/ATLAS` and OpenBlas is in :code:`/usr/local/openblas`. Matlab is installed in :code:`/usr/local/MATLAB/R2018b/`. In general (with MATLAB being an exception), each of these directories should contain a :code:`lib` and :code:`include` directory. You should point :code:`LDFLAGS` to the :code:`lib` directory and :code:`CPPFLAGS` to the :code:`include` directory. Thus, to configure with FFTW3 and OpenBlas, we run (as one line)

.. code-block:: bash

   CPPFLAGS="-I/usr/local/include -I/usr/local/openblas/include" \
                LDFLAGS="-L/usr/local/lib -L/usr/local/openblas/lib" \
   ./configure --with-fftw3 --with-fftw3-threads=omp 


Alternatively, you can export :code:`CPPFLAGS` and :code:`LDFLAGS`:
.. code-block:: bash

   export CPPFLAGS="-I/usr/local/include -I/usr/local/openblas/include"
   export LDFLAGS="-L/usr/local/lib -L/usr/local/openblas/lib"
   ./configure --with-fftw3 --with-fftw3-threads=omp 

mex bindings
==============

If we also we want to compile with the mex bindings, we need to add Matlab's :code:`lib` and :code:`include` directories. These are in non-stanard locations, so they must be added to :code:`CPPFLAGS` and :code:`LDFLAGS`. By default, the mex modules will get installed into :code:`${prefix}/lib/`, which is probably not what you want. Specify a different location with :code:`--with-mex-prefix=/path/to/mex`:

.. code-block:: bash

   export CPPFLAGS="-I/usr/local/include -I/usr/local/openblas/include \
                 -I/usr/local/MATLAB/R2018b/extern/include"
   export LDFLAGS="-L/usr/local/lib -L/usr/local/openblas/lib  \
                -L/usr/local/MATLAB/R2018b/bin/glnxa64"
   ./configure --with-fftw3 --with-fftw3-threads=omp --enable-mex \  
                --with-mex-prefix=/home/arnold/matlab/l1c

Note that on my system, the command :code:`mexext` is located in :code:`/usr/local/MATLAB/R2018b/bin/`, which is symlinked to :code:`/usr/local/bin/mexext`, which is on my path. If this is not the case, then in addition to above you can, e.g., :code:`export MEXEXT=mexa64`. You can get the appropriate value to export by typing :code:`mexext` at the matlab command prompt.


Python bindings
===============

To build the python bindings, use :code:`--enable-python`:

.. code-block:: bash

   export CPPFLAGS="-I/usr/local/include -I/usr/local/openblas/include"
   export LDFLAGS="-L/usr/local/lib -L/usr/local/openblas/lib"
   ./configure --with-fftw3 --with-fftw3-threads=omp --enable-python

Building python bindings is supported for Python 3 (tested with 3.5). The proper compilation and linking flags as well as the installation location are obtained from the python3 on your path (via distutils.sysconfig). On linux, the typical install location will default to something like :code:`/usr/lib/python3/dist-packages`. These values can be modified via the environmental variables:

.. code-block:: bash

   PYTHON_CPPFLAGS        # Should contain Python.h
   PYTHON_LIBS            # e.g., -lpython3.5m
   PYTHON_SITE_PKG_EXEC   # e.g., /home/user/.local/lib/python3.5/site-packages


Unit Tests
==============
Almost all of the test data is generated in python and saved as json files in :code:`$(build_dir)/test/test_data/`.
To run the test suite, execute 

.. code-block:: bash

   make check

For more information about the tests, see :ref:`testing`.


