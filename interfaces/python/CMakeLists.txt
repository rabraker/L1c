# Build python bindings

add_library(_l1cPy_module SHARED l1cPy.c)
set_target_properties(_l1cPy_module PROPERTIES PREFIX "")

target_link_libraries(_l1cPy_module
  PRIVATE l1c
  Python3::NumPy Python3::Module cjson)

list(APPEND PY_FILES
  __init__.py
  L1c.py
  test_L1cPy.py
  _examples/dct_example.py
  _examples/nesta_example.py
  _examples/l1c_pyplot_utils.py
  _examples/run_TV_denoise_bregman.py
  )

# Create the test data and setup python unit tests.

set(EXAMPLE_DATA_DEST ${CMAKE_CURRENT_BINARY_DIR}/example_img_data.json)
message("${EXAMPLE_DATA_DEST}")

foreach(py_file IN ITEMS ${PY_FILES})
  configure_file(${py_file} ${CMAKE_CURRENT_BINARY_DIR}/${py_file} COPYONLY)
endforeach()
set(L1C_PY_ROOT ${CMAKE_CURRENT_BINARY_DIR})
configure_file(_examples/l1c_py_init_path.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/_examples/l1c_py_init_path.py)

# This acts as a fixture to generate our test data. Below, we make the actual
# test depend on this target.
add_test(
  NAME setupExampleData
  COMMAND python3 ${CMAKE_SOURCE_DIR}/test/build_CS20NG_example_data.py ${EXAMPLE_DATA_DEST})

add_test(
  NAME test_python_interface
  COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/test_L1cPy.py)

set_tests_properties(setupExampleData PROPERTIES FIXTURES_SETUP example_data)
set_tests_properties(test_python_interface PROPERTIES FIXTURES_REQUIRED example_data)
