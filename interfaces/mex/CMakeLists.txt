
find_package(Octave REQUIRED COMPONENTS Development)

add_library(l1qc_dct MODULE
  l1c_mex_utils.c
  l1qc_dct_mex.c)
target_link_libraries(l1qc_dct
  Octave::Octave
  l1c::l1c)
set_target_properties(l1qc_dct PROPERTIES PREFIX ""
  SUFFIX ".mex")

add_library(nesta_dctTV MODULE
  l1c_mex_utils.c
  nesta_dctTV_mex.c)
target_link_libraries(nesta_dctTV
  Octave::Octave
  l1c::l1c)
set_target_properties(nesta_dctTV PROPERTIES PREFIX ""
  SUFFIX ".mex")

add_library(breg_anistropic_tv MODULE
  l1c_mex_utils.c
  breg_anistropic_TV_mex.c)
target_link_libraries(breg_anistropic_tv
  Octave::Octave
  l1c::l1c)
set_target_properties(breg_anistropic_tv PROPERTIES PREFIX ""
  SUFFIX ".mex")

# target_include_directories(l1qc_dct
#   PRIVATE "/usr/include/octave-4.2.2/octave")

set(_m_files
  get_data_dir.m
 breg_anistropic_TV.m
 l1qc_dct.m
 l1qc_dct_opts.m
 nesta_opts.m
 TestL1cMex.m
 _examples/cs20ng_grating.m
 _examples/l1c_mex_init_paths.m.in
 _examples/Makefile
 _examples/Makefile.am
 _examples/Makefile.in
 _examples/mu_path_mask.m
 _examples/nesta_example.m
 _examples/play_l1qc_lb_and_TV_example.m
 )

foreach(_file IN ITEMS ${_m_files})
  message(${_file})
  configure_file(${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file} COPYONLY)
endforeach()

add_subdirectory(test)
