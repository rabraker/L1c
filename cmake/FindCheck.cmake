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

set(CHECK_FOUND ${PC_Check_FOUND})
message("**** check found ${PC_Check_FOUND}")
message("**** check flags ${PC_Check_LDFLAGS}")
message("**** check flags ${PC_Check_LIBRARY_DIRS}")
message("**** check flags ${PC_Check_LINK_LIBRARIES}")

if( (Check_INCLUDE_DIR STREQUAL "check_INCLUDE_DIR-NOTOFUND"))
  set(CHECK_FOUND FALSE)
endif()


find_package_handle_standard_args(Check
  # DEFAULT_MSG
  REQUIRED_VARS Check_INCLUDE_DIR PC_Check_LDFLAGS
  VERSION_VAR PC_Check_VERSION)

mark_as_advanced(check_LIBRARY  check_INCLUDE_DIR)

if(cJSON_FOUND)
  add_library(Check::Check SHARED IMPORTED)

  set_target_properties(
    Check::Check
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PC_Check_INCLUDE_DIRS}"
    IMPORTED_LOCATION "${PC_Check_LINK_LIBRARIES}"

    )
  target_link_options(Check::Check INTERFACE
    ${PC_Check_LDFLAGS}
    )

endif()
