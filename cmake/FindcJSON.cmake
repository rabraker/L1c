# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
FindCJSON
---------

Finds the cJSON library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``cJSON_FOUND``
True if the system has the cJSON library.
``cJSON_VERSION``
The version of the cJSON library which was found.
``cJSON_INCLUDE_DIRS``
Include directories needed to use cJSON.
``cJSON_LIBRARIES``
Libraries needed to link to cJSON.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``cJSON_INCLUDE_DIR``
The directory containing ``foo.h``.
``cJSON_LIBRARY``
The path to the cJSON library.

#]=======================================================================]
include(FindPackageHandleStandardArgs)
find_package(PkgConfig)
pkg_check_modules(PC_cJSON libcjson)

set(cJSON_FOUND FALSE)

find_path(cJSON_INCLUDE_DIR
  NAMES cJSON.h
  PATHS ${PC_cJSON_INCLUDE_DIRS}
  PATH_SUFFIXES cjson
  )
find_library(cJSON_LIBRARY
  NAMES cjson
  PATHS ${PC_cJSON_LIBRARY_DIRS}
  )

find_package_handle_standard_args(cJSON DEFAULT_MSG
  cJSON_LIBRARY cJSON_INCLUDE_DIR)

mark_as_advanced(cJSON_LIBRARY  cJSON_INCLUDE_DIR)

if(cJSON_FOUND)
  add_library(cJSON::cJSON SHARED IMPORTED)

  set_target_properties(
    cJSON::cJSON
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PC_cJSON_INCLUDE_DIRS}"
    IMPORTED_LOCATION "${cJSON_LIBRARY}"
    )
endif()
