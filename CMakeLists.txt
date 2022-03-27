cmake_minimum_required(VERSION 3.15)

# Determine if am is built as a subproject (using add_subdirectory) or if it is
# the master project.
if(NOT DEFINED NDSU3LIB_MASTER_PROJECT)
  set(NDSU3LIB_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NDSU3LIB_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

project(ndsu3lib Fortran)

if(NDSU3LIB_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# ##############################################################################
# find external projects/dependencies
# ##############################################################################
set(SO3COEF_LIBRARY_OPTIONS "wigxjpf" "gsl")
set(SO3COEF_LIBRARY "gsl" CACHE STRING "SO(3) coefficient library to use")
set_property(CACHE SO3COEF_LIBRARY PROPERTY STRINGS ${SO3COEF_LIBRARY_OPTIONS})

if(SO3COEF_LIBRARY STREQUAL "wigxjpf")
  find_package(wigxjpf REQUIRED)
elseif(SO3COEF_LIBRARY STREQUAL "gsl")
  find_package(GSL REQUIRED)
else()
  message(FATAL_ERROR "Unknown SO(3) coefficient library ${SO3COEF_LIBRARY}")
endif()
message(STATUS "Found SO(3) coefficient library:  ${SO3COEF_LIBRARY}")


# ##############################################################################
# add sources
# ##############################################################################

set(ndsu3lib_sources
  ndsu3lib_recoupling.F90
  ndsu3lib_wigner_canonical.F90
  ndsu3lib_wigner_su3so3.F90
  c++wrappers.cpp
  c++wrappers.h
)

# construct lists of headers and sources
# define ndsu3lib library and add source files
add_library(ndsu3lib ${ndsu3lib_sources})
add_library(ndsu3lib::ndsu3lib ALIAS ndsu3lib)

# ##############################################################################
# configure options
# ##############################################################################

option(NDSU3LIB_OPENMP "Make ndsu3lib safe for use in OpenMP parallel regions" ON)

if(${NDSU3LIB_OPENMP})
  find_package(OpenMP REQUIRED)
  target_link_libraries(ndsu3lib OpenMP::OpenMP_Fortran)
endif()


## What is NDSU3LIB_QUAD?  Is it the intel compiler quad?
#
## define options 
set(NDSU3LIB_PRECISION_OPTIONS "double" "quad" "multi")
set(NDSU3LIB_PRECISION "double" CACHE STRING "NDSU3LIB numerical precision")
set_property(CACHE NDSU3LIB_PRECISION PROPERTY STRINGS ${NDSU3LIB_PRECISION_OPTIONS})

## set precision 
if(${NDSU3LIB_PRECISION} STREQUAL "double")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_DBL)
## [TODO: check correct] Assumes only use gsl quad precision if using gsl library for so3coefs 
elseif(${NDSU3LIB_PRECISION} STREQUAL "quad")
  if(SO3COEF_LIBRARY STREQUAL "gsl")
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD)
  endif()

## [TODO: check correct] Assumes only use gsl multi precision if using gsl library for so3coefs 
elseif(${NDSU3LIB_PRECISION} STREQUAL "multi")
  if(SO3COEF_LIBRARY STREQUAL "gsl")
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP)
  endif()
  
else()
  message(FATAL_ERROR "Precision option ${NDSU3LIB_PRECISION} not recognized")
endif()


# ##############################################################################
# define installation rules
# ##############################################################################

install(
  TARGETS ndsu3lib
  DESTINATION lib
  EXPORT ndsu3libTargets
)

install(
  EXPORT ndsu3libTargets
  NAMESPACE ndsu3lib::
  FILE ndsu3libTargets.cmake
  DESTINATION lib/cmake/ndsu3lib
)

include(CMakePackageConfigHelpers)
# generate the config file that is includes the exports
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ndsu3libConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/ndsu3lib"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ndsu3libConfigVersion.cmake"
  VERSION "${ndsu3lib_VERSION_MAJOR}.${ndsu3lib_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/ndsu3libConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/ndsu3libConfigVersion.cmake
        DESTINATION lib/cmake/ndsu3lib
)

export(
  EXPORT ndsu3libTargets
  NAMESPACE ndsu3lib::
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ndsu3libTargets.cmake"
)
