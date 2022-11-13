cmake_minimum_required(VERSION 3.15)

project(l1c
  VERSION "0.0.1"
  )
include(CheckCCompilerFlag)
include(CheckSymbolExists)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1)


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
list(APPEND L1C_CFLAGS_
  -O3
  -mtune=native
  -m64
  -march=native
  -fno-common
  -pedantic
  -Wall
  -Wextra
  -Wunused
  -Wcast-qual
  -Wcast-align
  -Winit-self
  -Wundef
  -Wshadow
  -Wstringop-truncation
  -Wswitch-enum
  -Wvla
  )
foreach(flag IN ITEMS ${L1C_CFLAGS_})
  string(REGEX REPLACE "[\-=]" "_" FLAG_VAR ${flag})
  check_c_compiler_flag(${flag} HAS${FLAG_VAR})
  if(${HAS${FLAG_VAR}})
    list(APPEND L1C_CFLAGS ${flag})
  endif()
endforeach()

option(USE_WERROR OFF)

if(USE_WERROR)
  add_compile_options(-Werror)
endif()

if(MSVC)
  add_compile_options(/W4 /Wx)
else()
  add_compile_options(
    ${L1C_CFLAGS}
    )
endif()

check_symbol_exists(_mm_malloc "xmmintrin.h" _HAVE_MM_MALLOC_)
check_symbol_exists(posix_memalign "stdlib.h" _HAVE_POSIX_MEMALIGN_)
if(_HAVE_POSIX_MEMALIGN_)
  # Unlocks posix_memalign in stdlib.h
  set(_POSIX_C_SOURCE 200112L)
endif()

include(${CMAKE_BINARY_DIR}/conan_paths.cmake)


set(BLAS_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(cJSON REQUIRED)

find_package(Python3
  REQUIRED
  COMPONENTS Interpreter NumPy Development)

add_subdirectory(src)
add_subdirectory(examples/c)
add_subdirectory(test)
add_subdirectory(interfaces/python)

option(ENABLE_MEX "Allow building the octave interface.")
if(ENABLE_MEX)
  add_subdirectory(interfaces/mex)
endif()


configure_file(
  "${PROJECT_SOURCE_DIR}/l1c_config.h.in"
  "${PROJECT_BINARY_DIR}/include/config.h"
  )
