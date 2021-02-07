# Build python bindings

add_library(_l1cPy_module SHARED l1cPy.c)
set_target_properties(_l1cPy_module PROPERTIES PREFIX "")

target_link_libraries(_l1cPy_module PUBLIC l1c Python3::Module cjson)

list(APPEND PY_FILES
  __init__.py
  L1c.py
  test_L1cPy.py
  _examples/dct_example.py
  # _examples/example_img_data.json
  # _examples/l1c_py_init_path.py
  _examples/l1c_pyplot_utils.py
  _examples/run_TV_denoise_bregman.py
  )

# Create the test data.
# $(test_data):$(abs_top_srcdir)/test/build_CS20NG_example_data.py
# python3 $< $@
set(EXAMPLE_DATA_DEST ${CMAKE_CURRENT_BINARY_DIR}/example_img_data.json)
message("${EXAMPLE_DATA_DEST}")

add_custom_command(
  OUTPUT ${EXAMPLE_DATA_DEST}
  COMMENT  "Generating example data for python bindings."
  VERBATIM
  COMMAND python3 ${CMAKE_SOURCE_DIR}/test/build_CS20NG_example_data.py ${EXAMPLE_DATA_DEST}
  )

add_custom_target(example_data ALL DEPENDS ${EXAMPLE_DATA_DEST})
# add_dependencies(_l1cPy_module example_data)

#

# add_custom_target(
#   example_data ALL
#   COMMAND python3 ${CMAKE_SOURCE_DIR}/test/build_CS20NG_example_data.py ${EXAMPLE_DATA_DEST}
#   # echo "*********** Echo Doing stuff"
#   DEPENDS ${CMAKE_SOURCE_DIR}/test/build_CS20NG_example_data.py
#   BYPRODUCTS ${EXAMPLE_DATA_DEST}
#   # DEPENDS ${EXAMPLE_DATA_DEST}
#   )


foreach(py_file IN ITEMS ${PY_FILES})
  # file(COPY ${py_file}
  #   DESTINATION ${CMAKE_BINARY_DIR}/py_file
  #   )
  # message("pyfile=${py_file}")
  configure_file(${py_file} ${CMAKE_CURRENT_BINARY_DIR}/${py_file} COPYONLY)
endforeach()
