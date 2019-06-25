.. _testing:

Testing Code
============


To run the all of the test suites, execute 

.. code-block:: bash

   make check

This will first build a set of test data and then execute multiple suites. The tests for :code:`libl1c` proper use :code:`libcheck`. Each of the bindings (python and matlab) have their own test suites. You can execute individual suites via

.. code-block:: bash

   make check TESTS='TEST_NAME'

where :code:`TEST_NAME=[test_matlab_mex.sh|test_python_interface.sh|test_l1c]`. Each of these test programs/scripts contains multiple suites which contain many cases. For each of these programs/scripts, you can run individual suites or cases:

.. code-block:: bash

   # Run a specific suite from the check suite
   CK_RUN_SUITE=nesta make check TESTS=test_l1c

   # Run a specific case from the check suite
   CK_RUN_SUITE=nesta make check TESTS=test_l1c

   #Run a specific case from the matlab tests
   MEX_RUN_CASE=check_nesta_dctTV make check TESTS='test_matlab_mex.sh'

   # Run a specific suite from the python tests
   PY_RUN_CASE=TestNesta make check TESTS='test_python_interface.sh'

   # Run a specific case from the python tests
   PY_RUN_CASE=TestNesta.test_nesta_errors make check TESTS='test_python_interface.sh'
   

By default, :code:`make check` will skip the memory leak test, which is very time consuming. To run this also, execute

.. code-block:: bash

   with_valgrind=yes make check


The mex test code uses a custom unit test framework, rather than the built-in matlab unittest framework. The goal here is to be able to use the same code to test an :code:`octave` interface, which as far as I can tell, is not compatible with the matlab unittest framework.
