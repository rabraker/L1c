#[=======================================================================[.rst:
FindFFTW3
-----------

This module finds the fftw3 library on the system.

Imported Targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` target:

``FFTW3::FFTW3``
  The fftw3 library, if found.

Result Variables
^^^^^^^^^^^^^^^^

The following variables are set:

``FFTW3_THREADS_ENABLED``
  If fftw3 was found with thread support.


#]=======================================================================]

include (CheckLibraryExists)
set(FFTW3_FOUND FALSE)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
set(CMAKE_REQUIRED_QUIET ${Threads_FIND_QUIETLY})

include (CheckIncludeFile)
include (CheckCSourceCompiles)


find_package(PkgConfig)
pkg_check_modules(PC_FFTW3 fftw3)

find_path(FFTW3_INCLUDE_DIRS
  NAMES fftw3.h
  PATHS ${PC_FFTW3_INCLUDE_DIRS}
  PATH_SUFFIXES Foo
  )

find_library(FFTW3_LIBRARIES
  NAMES fftw3_omp fftw3_threads fftw3
  PATHS ${PC_FFTW3_LIBRARY_DIRS}
  PATH_SUFFIXES Foo
  )



# simple fftw3 thread test code
set(FFTW3_C_TEST_SOURCE [====[
#include <fftw3.h>

/* Override any GCC internal prototype to avoid an error.
   Use char because int might match the return type of a GCC
   builtin and then its argument prototype would still apply.  */
#ifdef __cplusplus
extern "C"
#endif

int main ()
{
return fftw_init_threads ();
  ;
  return 0;
}

]====])


find_package_handle_standard_args(FFTW3 DEFAULT_MSG
  FFTW3_INCLUDE_DIRS FFTW3_LIBRARIES)

# mark_as_advanced(FFTW3_INCLUDE_DIRS  FFTW3_LIBRARIES)

unset(FFTW3_INCLUDE_DIRS)
unset(FFTW3_LIBRARIES)

# Check if we have threads.
set(CMAKE_REQUIRED_LIBRARIES "${FFTW3_LIBRARIES}")
CHECK_C_SOURCE_COMPILES(
  "${FFTW3_C_TEST_SOURCE}"
  CMAKE_HAVE_LIB_FFTW3_THREADS
  )

add_library(FFTW3::FFTW3 SHARED IMPORTED)
set_target_properties(
  FFTW3::FFTW3
  PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}"
  IMPORTED_LOCATION "${FFTW3_LIBRARIES}"
  )

if(CMAKE_HAVE_LIB_FFTW3_THREADS)
  set(HAVE_FFTW3_THREADS TRUE)
else()
  unset(HAVE_FFTW3_THREADS)
endif()
