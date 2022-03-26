
add_executable(l1qc_dct_c l1qc_dct_c.c)
target_link_libraries(l1qc_dct_c
  PRIVATE l1c::l1c cJSON::cJSON

  )
