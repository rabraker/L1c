
set(L1C_SOURCES
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
  vcl/vcl_math.cpp
  )

add_library(l1c SHARED "${L1C_SOURCES}")
find_package(OpenMP)

target_link_libraries(l1c
  PRIVATE
  ${L1C_PUBLIC_LIBS}
  FFTW3::FFTW3
  OpenMP::OpenMP_C
  )

target_include_directories(l1c PUBLIC include PRIVATE include)
