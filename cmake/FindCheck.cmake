# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindCJSON
---------

Finds the check unittest library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``check_FOUND``
True if the system has the check library.
``check_VERSION``
The version of the cJSON library which was found.

Cache Variables
^^^^^^^^^^^^^^^


#]=======================================================================]
set(check_FOUND FALSE)

include(FindPackageHandleStandardArgs)
find_package(PkgConfig)
pkg_check_modules(PC_Check check)

set(Check_VERSION ${PC_Check_VERSION})

find_path(Check_INCLUDE_DIR
  NAMES check.h
  PATHS ${PC_Check_INCLUDE_DIRS}
  )

# PC_Check_LINK_LIBRARIES can contain a comma seperated list of all the other
# transitive link libraries, like /usr/lib/x86_64-linux-gnu/librt.so. Look for
# the one item with libcheck in its name to use for the IMPORTED_LOCATION.
foreach(lib IN ITEMS ${PC_Check_LINK_LIBRARIES})
  STRING(FIND ${lib} libcheck _FOUND)
  if(NOT ("${_FOUND}" STREQUAL "-1"))
    set(Check_LOCATION ${lib})
    break()
  endif()

endforeach()


find_package_handle_standard_args(Check
  REQUIRED_VARS Check_INCLUDE_DIR PC_Check_LDFLAGS Check_LOCATION
  VERSION_VAR PC_Check_VERSION)

mark_as_advanced(check_LIBRARY  check_INCLUDE_DIR)

if(Check_FOUND)
  add_library(Check::Check SHARED IMPORTED)

  set_target_properties(
    Check::Check
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PC_Check_INCLUDE_DIRS}"
    IMPORTED_LOCATION ${Check_LOCATION}
    )

  target_link_libraries(Check::Check INTERFACE
    ${PC_Check_LDFLAGS}
    )

endif()
