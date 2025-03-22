# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
Findwigxjpf
-----------

Finds the wigxjpf library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``wigxjpf::wigxjpf``
  The wigxjpf library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``wigxjpf_FOUND``
  True if the system has the wigxjpf library.
``wigxjpf_INCLUDE_DIRS``
  Include directories needed to use wigxjpf.
``wigxjpf_LIBRARIES``
  Libraries needed to link to wigxjpf.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``wigxjpf_INCLUDE_DIR``
  The directory containing ``fwigxjpf.mod``.
``wigxjpf_LIBRARY``
  The path to the wigxjpf library.

#]=======================================================================]

find_path(wigxjpf_INCLUDE_DIR
  NAMES fwigxjpf.mod
  PATHS ${WIGXJPF_INC_DIR}
  PATH_SUFFIXES mod wigxjpf
)

find_library(wigxjpf_LIBRARY
  NAMES wigxjpf libwigxjpf libwigxjpf.a wigxjpf.a
  PATHS ${WIGXJPF_LIB_DIR}
  PATH_SUFFIXES lib wigxjpf
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(wigxjpf
  REQUIRED_VARS
    wigxjpf_LIBRARY
    wigxjpf_INCLUDE_DIR
)

if(wigxjpf_FOUND)
  set(wigxjpf_LIBRARIES ${wigxjpf_LIBRARY})
  set(wigxjpf_INCLUDE_DIRS ${wigxjpf_INCLUDE_DIR})
  set(wigxjpf_DEFINITIONS )
endif()

if(wigxjpf_FOUND AND NOT TARGET wigxjpf::wigxjpf)
  add_library(wigxjpf::wigxjpf STATIC IMPORTED)
  set_target_properties(wigxjpf::wigxjpf PROPERTIES
    IMPORTED_LOCATION "${wigxjpf_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${wigxjpf_INCLUDE_DIR}"
  )
endif()

mark_as_advanced(
  wigxjpf_INCLUDE_DIR
  wigxjpf_LIBRARY
)
