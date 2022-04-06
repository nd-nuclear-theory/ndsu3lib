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

project(ndsu3lib LANGUAGES Fortran CXX)

if(NDSU3LIB_MASTER_PROJECT AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


# define module directory
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mod)
###############################################################################
# Find external projects/dependencies
###############################################################################
set(SO3COEF_LIBRARY_OPTIONS "wigxjpf" "gsl")
set(SO3COEF_LIBRARY "gsl" CACHE STRING "SO(3) coefficient library to use")
set_property(CACHE SO3COEF_LIBRARY PROPERTY STRINGS ${SO3COEF_LIBRARY_OPTIONS})

if(NDSU3LIB_MASTER_PROJECT)
  include(retrieve.cmake)
endif()

if(SO3COEF_LIBRARY STREQUAL "wigxjpf")
  retrieve_package(wigxjpf::wigxjpf wigxjpf)

elseif(SO3COEF_LIBRARY STREQUAL "gsl")
  find_package(GSL REQUIRED)

else()
  message(FATAL_ERROR "Unknown SO(3) coefficient library ${SO3COEF_LIBRARY}")

endif()

message(STATUS "Using SO(3) coefficient library:  ${SO3COEF_LIBRARY}")

###############################################################################
# Set precision
###############################################################################
set(NDSU3LIB_PRECISION_OPTIONS "double" "quad" "multi")
set(NDSU3LIB_PRECISION "double" CACHE STRING "NDSU3LIB numerical precision")
set_property(CACHE NDSU3LIB_PRECISION PROPERTY STRINGS ${NDSU3LIB_PRECISION_OPTIONS})

## mpfun20-fortran
if(${NDSU3LIB_PRECISION} STREQUAL "multi")
  retrieve_package(mpfun20:mpfun20 mpfun20)
endif()

###############################################################################
# Add library
###############################################################################
set(ndsu3lib_sources_F
  ndsu3lib_recoupling.F90
  ndsu3lib_wigner_canonical.F90
  ndsu3lib_wigner_su3so3.F90
)
set(ndsu3lib_sources_HEADERS  c++wrappers.h)
set(ndsu3lib_sources_CPP c++wrappers.cpp)
set(ndsu3lib_sources ${ndsu3lib_sources_F} ${ndsu3lib_sources_CPP})

add_library(ndsu3lib ${ndsu3lib_sources})
add_library(ndsu3lib::ndsu3lib ALIAS ndsu3lib)

target_sources(ndsu3lib INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/c++wrappers.h>)

## MKL
find_package(BLAS)
target_link_libraries(ndsu3lib PUBLIC BLAS::BLAS)


###############################################################################
#configure options
###############################################################################
if(${SO3COEF_LIBRARY} STREQUAL "gsl")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_RACAH_GSL)
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_WSO3_GSL)
  target_link_libraries(ndsu3lib PUBLIC GSL::gsl)

elseif(${SO3COEF_LIBRARY} STREQUAL "wigxjpf")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_RACAH_WIGXJPF)
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_WSO3_WIGXJPF)
  target_link_libraries(ndsu3lib PUBLIC wigxjpf::wigxjpf)

endif()

## set precision 
if(${NDSU3LIB_PRECISION} STREQUAL "double")
  target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_DBL)

elseif(${NDSU3LIB_PRECISION} STREQUAL "quad")
  
  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_QUAD)
  endif()

elseif(${NDSU3LIB_PRECISION} STREQUAL "multi")
  target_link_libraries(ndsu3lib PUBLIC mpfun20::mpfun20)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP_GNU)
  else()
    target_compile_definitions(ndsu3lib PRIVATE NDSU3LIB_MP)
  endif()
  
else()
  message(FATAL_ERROR "Precision option ${NDSU3LIB_PRECISION} not recognized")

endif()

message(STATUS "Precision: ${NDSU3LIB_PRECISION}")


## OpenMP Not currently implemented? 
option(NDSU3LIB_OPENMP "Make ndsu3lib safe for use in OpenMP parallel regions" ON)

if(${NDSU3LIB_OPENMP})
  find_package(OpenMP REQUIRED)
  target_link_libraries(ndsu3lib PUBLIC OpenMP::OpenMP_Fortran)
endif()


###############################################################################
#define include directory
###############################################################################
target_include_directories(
  ndsu3lib
  PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/mod>
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
file(GLOB mod_list ${CMAKE_Fortran_MODULE_DIRECTORY}/ndsu3lib*.mod)
foreach(mod IN LISTS mod_list)
  install(FILES ${CMAKE_Fortran_MODULE_DIRECTORY}/mod DESTINATION include/ndsu3lib)
endforeach()


install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/c++wrappers.h DESTINATION include/ndsu3lib)

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

add_executable(ndsu3lib_example_cpp ndsu3lib_example.cpp)
target_link_libraries(ndsu3lib_example_cpp ndsu3lib::ndsu3lib)
set_target_properties(ndsu3lib_example_cpp PROPERTIES LINKER_LANGUAGE CXX)


if(NDSU3LIB_MASTER_PROJECT)
  add_custom_target(tests)
  add_dependencies(tests ndsu3lib_example)
  add_dependencies(tests ndsu3lib_example_cpp)
endif()

