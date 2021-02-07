
add_executable(l1qc_dct l1qc_dct_c.c)
target_link_libraries(l1qc_dct
  PRIVATE l1c FFTW3::FFTW3 ${L1C_PUBLIC_LIBS}
  )
