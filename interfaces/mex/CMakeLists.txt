find_package(Octave REQUIRED COMPONENTS Development)

add_library(l1qc_dct MODULE l1c_mex_utils.c l1qc_dct_mex.c)
target_link_libraries(l1qc_dct Octave::Octave l1c::l1c)
set_target_properties(l1qc_dct PROPERTIES PREFIX "" SUFFIX ".mex")

add_library(nesta_dctTV MODULE l1c_mex_utils.c nesta_dctTV_mex.c)
target_link_libraries(nesta_dctTV Octave::Octave l1c::l1c)
set_target_properties(nesta_dctTV PROPERTIES PREFIX "" SUFFIX ".mex")

add_library(breg_anistropic_tv MODULE l1c_mex_utils.c breg_anistropic_TV_mex.c)
target_link_libraries(breg_anistropic_tv Octave::Octave l1c::l1c)
set_target_properties(breg_anistropic_tv PROPERTIES PREFIX "" SUFFIX ".mex")

set(_m_files
    breg_anistropic_TV.m
    l1qc_dct.m
    l1qc_dct_opts.m
    nesta_opts.m
    examples/cs20ng_grating.m
    examples/l1c_mex_init_paths.m.in
    examples/mu_path_mask.m
    examples/nesta_example.m
    examples/play_l1qc_lb_and_TV_example.m)

foreach(_file IN ITEMS ${_m_files})
  message(${_file})
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

set(_mex_test_srcs
    TMU_mex_assert_2Darray_with_size.c
    TMU_mex_assert_num_inputs.c
    TMU_mex_assert_num_outputs.c
    TMU_mex_assert_scalar_struct.c
    TMU_mex_get_double_array_or_fail.c
    TMU_mex_get_double_from_struct_or_fail.c
    TMU_mex_get_double_scalar_or_fail.c
    TMU_mex_get_int_array_or_fail.c)

foreach(_test_src IN ITEMS ${_mex_test_srcs})
  string(REGEX REPLACE "\.c$" "" _test_src_name ${_test_src})
  message(${_test_src_name})

  add_library(${_test_src_name} MODULE test/${_test_src} l1c_mex_utils.c)
  target_link_libraries(${_test_src_name} Octave::Octave l1c::l1c)
  set_target_properties(${_test_src_name} PROPERTIES PREFIX "" SUFFIX ".mex")
  target_include_directories(${_test_src_name}
                             PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/)

endforeach()

set(_m_test_files test/L1cMexTesting.m test/test_mex_interface.m
                  test/test_mex_utils.m test/clrs.m)

foreach(_file IN ITEMS ${_m_test_files})
  message(${_file})
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

add_test(NAME mex_nesta
         COMMAND octave --path=${CMAKE_CURRENT_BINARY_DIR}/examples
                 ${CMAKE_CURRENT_BINARY_DIR}/examples/nesta_example.m)
add_test(
  NAME play_l1qc_lb_and_TV_example
  COMMAND octave --path=${CMAKE_CURRENT_BINARY_DIR}/examples
          ${CMAKE_CURRENT_BINARY_DIR}/examples/play_l1qc_lb_and_TV_example.m)

add_test(NAME test_mex_utils
         COMMAND octave --path=${CMAKE_CURRENT_BINARY_DIR}/test
                 ${CMAKE_CURRENT_BINARY_DIR}/test/test_mex_utils.m)

add_test(NAME test_mex_interface
         COMMAND octave --path=${CMAKE_CURRENT_BINARY_DIR}/test
                 ${CMAKE_CURRENT_BINARY_DIR}/test/test_mex_interface.m)
