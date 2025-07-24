cmake_minimum_required(VERSION 3.15)

# Determine if am is built as a subproject (using add_subdirectory)
#  or if it is the master project.
if(NOT DEFINED NDSU3LIB_MASTER_PROJECT)
  set(NDSU3LIB_MASTER_PROJECT OFF)
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(NDSU3LIB_MASTER_PROJECT ON)
    message(STATUS "CMake version: ${CMAKE_VERSION}")
  endif()
endif()

project(ndsu3lib LANGUAGES Fortran C CXX)

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 11)

if(NDSU3LIB_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


# define module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mod)
###############################################################################
# Find external projects/dependencies
###############################################################################
set(SU2COEF_LIBRARY_OPTIONS "wigxjpf" "gsl")
set(SU2COEF_LIBRARY "gsl" CACHE STRING "SO(3) coefficient library to use")
set_property(CACHE SU2COEF_LIBRARY PROPERTY STRINGS ${SU2COEF_LIBRARY_OPTIONS})


include(fetch_declarations.cmake)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR})

option(FETCH "Enable fetching of WIGXJPF or MPFUN20 from github.com/nd-nuclear-theory")

if(SU2COEF_LIBRARY STREQUAL "wigxjpf")
  message(STATUS "Building ndsu3lib with wigxjpf support")

  if(NOT TARGET wigxjpf::wigxjpf)
    find_package(wigxjpf)
  endif()
  ## If still not found, fetch
  if(NOT TARGET wigxjpf::wigxjpf AND FETCH)
    message(STATUS "Package wigxjpf not found. Fetching...")
    FetchContent_MakeAvailable(wigxjpf)
  endif()
  if(NOT TARGET wigxjpf::wigxjpf)
    message(FATAL_ERROR "WIGXJPF not found nor fetched")
  endif()

elseif(SU2COEF_LIBRARY STREQUAL "gsl")
  message(STATUS "Building ndsu3lib with gsl support")
  find_package(GSL REQUIRED)

else()
  message(FATAL_ERROR "Unknown SO(3) coefficient library ${SU2COEF_LIBRARY}")

endif()

message(STATUS "Using SO(3) coefficient library:  ${SU2COEF_LIBRARY}")

###############################################################################
# Set precision
###############################################################################
set(PRECISION_OPTIONS "double" "quad" "multi" "multiquad")
set(PRECISION "double" CACHE STRING "NDSU3LIB numerical precision")
set_property(CACHE PRECISION PROPERTY STRINGS ${PRECISION_OPTIONS})

## mpfun20-fortran
if(${PRECISION} STREQUAL "multi" OR ${PRECISION} STREQUAL "multiquad")
  message(STATUS "Building ndsu3lib with mpfun2020 support")
  if(NOT TARGET mpfun20)
    find_package(mpfun20)
  endif()
  ## If still not found, fetch
  if(NOT TARGET mpfun20 AND FETCH)
    message(STATUS "Package mpfun20 not found. Fetching...")
    FetchContent_MakeAvailable(mpfun20)
  endif()
  if(NOT TARGET mpfun20)
    message(FATAL_ERROR "MPFUN20 not found nor fetched")
  endif()
endif()

###############################################################################
# Add library
###############################################################################
set(ndsu3lib_sources_F
  ndsu3lib_coupling_canonical.F90
  ndsu3lib_coupling_su3so3.F90
  ndsu3lib_recoupling.F90
  ndsu3lib_tools.F90
)
set(ndsu3lib_sources_HEADERS  ndsu3lib.h)
set(ndsu3lib_sources ${ndsu3lib_sources_F})

add_library(ndsu3lib ${ndsu3lib_sources})
add_library(ndsu3lib::ndsu3lib ALIAS ndsu3lib)

target_sources(ndsu3lib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ndsu3lib.h>)

## Lapack
#find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#target_link_libraries(ndsu3lib PUBLIC BLAS::BLAS)
target_link_libraries(ndsu3lib PUBLIC LAPACK::LAPACK)


###############################################################################
#configure options
###############################################################################
if(${SU2COEF_LIBRARY} STREQUAL "gsl")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_RACAH_GSL)
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_WSO3_GSL)
  target_link_libraries(ndsu3lib PUBLIC GSL::gsl)

elseif(${SU2COEF_LIBRARY} STREQUAL "wigxjpf")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_RACAH_WIGXJPF)
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_WSO3_WIGXJPF)
  target_link_libraries(ndsu3lib PUBLIC wigxjpf::wigxjpf)

endif()

## set precision
if(${PRECISION} STREQUAL "double")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_DBL)

elseif(${PRECISION} STREQUAL "quad")

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU OR CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD)
  endif()

elseif(${PRECISION} STREQUAL "multi")
  target_link_libraries(ndsu3lib PUBLIC mpfun20)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU OR CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP)
  endif()

elseif(${PRECISION} STREQUAL "multiquad")
  target_link_libraries(ndsu3lib PUBLIC mpfun20)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU OR CMAKE_Fortran_COMPILER_ID STREQUAL Flang)
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD_GNU NDSU3LIB_MP_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD NDSU3LIB_MP)
  endif()

else()
  message(FATAL_ERROR "Precision option ${PRECISION} not recognized")

endif()


message(STATUS "Precision: ${PRECISION}")


option(OPENMP "Make ndsu3lib safe for use in OpenMP parallel regions")

## this is needed for Cray compiler
if(CMAKE_Fortran_COMPILER_ID STREQUAL Cray)
  if(${OPENMP})
    add_link_options(-Wl,--allow-multiple-definition -fopenmp)
  else()
    add_link_options(-Wl,--allow-multiple-definition)
  endif()
endif()

if(${OPENMP})
  find_package(OpenMP REQUIRED)
  target_link_libraries(ndsu3lib PUBLIC OpenMP::OpenMP_Fortran)
endif()


###############################################################################
#define include directory
###############################################################################
target_include_directories(
  ndsu3lib
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  $<INSTALL_INTERFACE:include/ndsu3lib>
)

###############################################################################
#define installation rules
###############################################################################
install(
  TARGETS ndsu3lib
  DESTINATION lib
  EXPORT ndsu3libTargets
)

# Copy .mod files to include directory
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY} DESTINATION include/ndsu3lib FILES_MATCHING PATTERN "*.mod")
# Copy header to include directory
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/ndsu3lib.h DESTINATION include/ndsu3lib)

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

###############################################################################
#define tests
###############################################################################
add_executable(ndsu3lib_example ndsu3lib_example.f90)
target_link_libraries(ndsu3lib_example ndsu3lib::ndsu3lib)
set_target_properties(ndsu3lib_example PROPERTIES LINKER_LANGUAGE Fortran)

add_executable(ndsu3lib_example_cpp ndsu3lib_example_cpp.cpp)
target_link_libraries(ndsu3lib_example_cpp ndsu3lib::ndsu3lib)
set_target_properties(ndsu3lib_example_cpp PROPERTIES LINKER_LANGUAGE CXX)

add_executable(ndsu3lib_example_c ndsu3lib_example_c.c)
target_link_libraries(ndsu3lib_example_c ndsu3lib::ndsu3lib)
set_target_properties(ndsu3lib_example_c PROPERTIES LINKER_LANGUAGE C)

if(NDSU3LIB_MASTER_PROJECT)
  add_custom_target(tests)
  add_dependencies(tests ndsu3lib_example)
  add_dependencies(tests ndsu3lib_example_cpp)
  add_dependencies(tests ndsu3lib_example_c)
endif()

