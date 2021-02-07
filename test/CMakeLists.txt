#
# Tests
#

find_package(Check REQUIRED)

add_executable(test_l1c.test
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

target_link_libraries(test_l1c.test
  PRIVATE
  l1c::l1c
  m
  Check::Check
  BLAS::BLAS
  cJSON::cJSON
)

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

set(_test_data_dir ${CMAKE_CURRENT_BINARY_DIR}/test_data)

# Before the scripts can run, create the test/test_data folder.
add_test(
  NAME setupDataDir
  COMMAND mkdir -p ${_test_data_dir})


foreach(test_file IN ITEMS ${TEST_DATA_SCRIPT_SOURCES})
  set(test_file_ ${CMAKE_CURRENT_SOURCE_DIR}/${test_file})
  get_filename_component(_name ${test_file_} NAME)
  add_test(
    NAME ${_name}
    COMMAND python3 ${test_file_})

  set_tests_properties(${_name} PROPERTIES FIXTURES_SETUP ${_name}_data)
  set_tests_properties(${test_prog} PROPERTIES ENVIRONMENT
    "TEST_DATA_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_data;ABS_TOP_BUILDDIR=${CMAKE_BINARY_DIR}")
endforeach()


add_test(
  NAME test_l1c.test
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/test_l1c.test)

set_tests_properties(test_l1c.test PROPERTIES ENVIRONMENT
  "TEST_DATA_DIR=${CMAKE_CURRENT_BINARY_DIR}/test_data;ABS_TOP_BUILDDIR=${CMAKE_BINARY_DIR}")

foreach(test_file IN ITEMS ${TEST_DATA_SCRIPT_SOURCES})
  set_tests_properties(test_l1c.test
    PROPERTIES FIXTURES_REQUIRED ${test_file}_data)
endforeach()

set(_verbose 1)
set(_data_path ${_test_data_dir}/example_img_data127.json)
add_test(
  NAME test_l1qc_dct_c_example
  COMMAND ${CMAKE_BINARY_DIR}/examples/c/l1qc_dct_c ${_data_path} ${_verbose})
