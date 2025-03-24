# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file LICENSE.rst or https://cmake.org/licensing for details.

#[=======================================================================[.rst:
Findmpfun20
-----------

Finds the mpfun20 library.

Imported Targets
^^^^^^^^^^^^^^^^

This module provides the following imported targets, if found:

``mpfun20``
  The mpfun20 library

Result Variables
^^^^^^^^^^^^^^^^

This will define the following variables:

``mpfun20_FOUND``
  True if the system has the mpfun20 library.
``mpfun20_INCLUDE_DIRS``
  Include directories needed to use mpfun20.

Cache Variables
^^^^^^^^^^^^^^^

The following cache variables may also be set:

``mpfun20_INCLUDE_DIR``
  The directory containing mpfun20 *.mod files.

#]=======================================================================]

find_path(mpfun20_INCLUDE_DIR
  NAMES mpmodule.o mpfuna.o mpfunb.o mpfunc.o mpfund.o mpfune.o mpfunf.o mpfung2.o mpfunh2.o mpmask13.o secondmpmodule.mod mpfuna.mod mpfunb.mod mpfunc.mod mpfund.mod mpfune.mod mpfunf.mod mpfung.mod mpfunh.mod
  PATHS ${MPFUN20_DIR}
  PATH_SUFFIXES mod mpfun20 inc include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(mpfun20
  REQUIRED_VARS
    mpfun20_INCLUDE_DIR
)

if(mpfun20_FOUND)
  set(mpfun20_INCLUDE_DIRS ${mpfun20_INCLUDE_DIR})
  set(mpfun20_DEFINITIONS )
endif()

if(mpfun20_FOUND AND NOT TARGET mpfun20)
  add_library(mpfun20
    ${mpfun20_INCLUDE_DIR}/mpmodule.o
    ${mpfun20_INCLUDE_DIR}/mpfuna.o
    ${mpfun20_INCLUDE_DIR}/mpfunb.o
    ${mpfun20_INCLUDE_DIR}/mpfunc.o
    ${mpfun20_INCLUDE_DIR}/mpfund.o
    ${mpfun20_INCLUDE_DIR}/mpfune.o
    ${mpfun20_INCLUDE_DIR}/mpfunf.o
    ${mpfun20_INCLUDE_DIR}/mpfung2.o
    ${mpfun20_INCLUDE_DIR}/mpfunh2.o
    ${mpfun20_INCLUDE_DIR}/mpmask13.o
    ${mpfun20_INCLUDE_DIR}/second.o
  )
#  add_library(mpfun20::mpfun20 OBJECT IMPORTED)
  set_target_properties(mpfun20 PROPERTIES
    LINKER_LANGUAGE Fortran
#    IMPORTED_LOCATION "${mpfun20_INCLUDE_DIR}"
    INTERFACE_INCLUDE_DIRECTORIES "${mpfun20_INCLUDE_DIR}"
  )
  install(
    TARGETS mpfun20
    DESTINATION lib
    EXPORT mpfun20Targets
  )
  install(
    EXPORT mpfun20Targets
    NAMESPACE mpfun20::
    FILE mpfun20Targets.cmake
    DESTINATION lib/cmake/mpfun20
  )
  export(
    EXPORT mpfun20Targets
    NAMESPACE mpfun20::
    FILE "${CMAKE_CURRENT_BINARY_DIR}/mpfun20Targets.cmake"
  )
endif()

mark_as_advanced(
  mpfun20_INCLUDE_DIR
)
