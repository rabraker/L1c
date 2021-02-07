#
# Tests


list(APPEND TEST_L1C_TEST_SOURCES
  TEST_cgsolve.c
  TEST_l1c_math.c
  TEST_l1c_memory.c
  TEST_l1qc_newton.c
  TEST_vcl_math.c
  TEST_bregman.c
  TEST_TV.c
  TEST_linesearch.c
  TEST_l1c_logging.c
  TEST_nesta.c
  TEST_matrix_transforms.c
  test_l1c.c
  TEST_dct.c
  TEST_dct2.c
  TEST_dct_tv.c
  )
set(test_l1c_test_INCLUDES
  $(top_srcdir)/include/bregman.h
  $(top_srcdir)/include/check_utils.h
  $(top_srcdir)/include/json_utils.h
  $(top_srcdir)/include/l1c.h
  $(top_srcdir)/include/l1c_math.h
  $(top_srcdir)/include/l1c_timing.h
  $(top_srcdir)/include/l1qc_newton.h
  $(top_srcdir)/include/linesearch.h
  $(top_srcdir)/include/nesta.h
  $(top_srcdir)/include/test_constants.h
  $(top_srcdir)/include/TEST.h
  $(top_srcdir)/include/TV.h
  $(top_srcdir)/include/vcl_math.h
  )


add_executable(test_l1c.test
  ${TEST_L1C_TEST_SOURCES})

target_link_libraries(test_l1c.test PUBLIC check ${L1C_PUBLIC_LIBS})
target_link_libraries(test_l1c.test PRIVATE l1c)

# We have scripts that generate test data, and scripts that run tests.
# All of them get copied into the build ./test folder.

#========================================================== #
# Test data generating scripts
list(APPEND TEST_DATA_SCRIPT_SOURCES
  build_CS20NG_example_data.py
  build_dct2_data.py
  build_dct_data.py
  build_dctTV_data.py
  build_cgsolve_data.py
  build_l1qc_data.py
  build_nesta_data.py
  build_TV_data.py
  build_bregman_data.py
  build_matrix_transform_data.py
  )

list(APPEND TEST_DATA_PY_SOURCES
  ${TEST_DATA_SCRIPT_SOURCES}
  L1cTestDataUtils.py
  test_runner_utils.sh
  )

# Copy all the data generating scripts into the build tree.
foreach(test_file IN ITEMS ${TEST_DATA_PY_SOURCES})
  configure_file(${test_file} ${CMAKE_CURRENT_BINARY_DIR}/${test_file} COPYONLY)
  list(APPEND TEST_DATA_SCRIPTS ${CMAKE_CURRENT_BINARY_DIR}/${test_file})
endforeach()

# Before the scripts can run, create the test/test_data folder.
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/test_data
  COMMENT "Creating test_data dir"
  VERBATIM
  COMMAND mkdir ${CMAKE_CURRENT_BINARY_DIR}/test_data
  )

foreach(test_file IN ITEMS ${TEST_DATA_SCRIPT_SOURCES})
  set(test_file_ ${CMAKE_CURRENT_BINARY_DIR}/${test_file})
  get_filename_component(_name ${test_file_} NAME)
  add_custom_command(
    OUTPUT "${test_file}.log"
    COMMENT "Generating test data: ${test_file_}"
    VERBATIM
    COMMAND python3 ${test_file_} > "${test_file_}.log"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/test_data
    )
  list(APPEND TEST_DATA_LOGS "${test_file_}.log")
endforeach()

# You can run `make test_data` to run the scripts.
add_custom_target(
  test_data
  DEPENDS ${TEST_DATA_LOGS}
  )

#========================================================== #
# Test running scripts
list(APPEND TEST_PROG_SOURCES
  test_l1qc_dct_c.sh
  valgrind_wrapper.sh
  )

if(ENABLE_MATLAB)
  list(APPEND TEST_PROG_SOURCES test_matlab_mex.sh)
endif()

if(ENABLE_PYTHON)
  list(APPEND TEST_PROG_SOURCES test_python_interface.sh)
endif()

# Copy all the test scripts into the build tree.
foreach(test_file IN ITEMS ${TEST_PROG_SOURCES})
  configure_file(${test_file} ${CMAKE_CURRENT_BINARY_DIR}/${test_file} COPYONLY)
  list(APPEND TEST_PROGS ${CMAKE_CURRENT_BINARY_DIR}/${test_file})
endforeach()

list(APPEND TEST_PROGS ${CMAKE_CURRENT_BINARY_DIR}/test_l1c.test)

# add_dependencies(_l1cPy_module example_data)

foreach(test_prog IN ITEMS ${TEST_PROGS})
  # message(${test_prog})
  set(test_log "${test_prog}.log")
  add_custom_command(
    COMMENT "checking ${test_prog}"
    COMMAND ${CMAKE_COMMAND} -E env TEST_DATA_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_data
    WITH_PYTHON=${ENABLE_PYTHON} WITH_MEX=${ENABLE_MATLAB}
    ABS_TOP_BUILDDIR=${CMAKE_BINARY_DIR}
    "${CMAKE_SOURCE_DIR}/test-driver"
         --log-file ${test_prog}.log
         --trs-file ${test_prog}.trs
         --color-tests yes
         --test-name ${test_prog} ${test_prog}
      OUTPUT "${test_log}.none"
    )
  list(APPEND TESTS ${test_log}.none)
endforeach()

add_custom_target(
  test
  DEPENDS ${TESTS} test_data
  )
