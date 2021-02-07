find_package(OpenMP)

add_library(l1c SHARED
  l1qc_newton.c
  l1c_math.c
  cgsolve.c
  json_utils.c
  TV.c
  bregman.c
  l1qc_dct.c
  linesearch.c
  l1c_logging.c
  l1c_memory.c
  matrix_transforms.c
  dct_transforms.c
  dct1.c
  dct2.c
  dct_tv_transforms.c
  nesta.c
  vcl/vcl_math.cpp)

add_library(l1c::l1c ALIAS l1c)

target_link_libraries(l1c
  PRIVATE
  BLAS::BLAS
  cJSON::cJSON
  FFTW3::FFTW3
  OpenMP::OpenMP_C)

target_include_directories(l1c
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${PROJECT_BINARY_DIR}/include)
