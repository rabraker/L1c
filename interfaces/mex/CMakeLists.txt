
find_package(Octave REQUIRED COMPONENTS Development)

add_library(l1qc_dct MODULE
  l1c_mex_utils.c
  l1qc_dct_mex.c)
target_link_libraries(l1qc_dct
  Octave::Octave
  l1c::l1c)
set_target_properties(l1qc_dct PROPERTIES PREFIX ""
  SUFFIX ".mex")
# target_include_directories(l1qc_dct
#   PRIVATE "/usr/include/octave-4.2.2/octave")
