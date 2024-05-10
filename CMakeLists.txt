# The whole cmake structure of this project is
# heavily borrowing from
# https://github.com/dftd4/dftd4/

cmake_minimum_required(VERSION 3.18)

# Load modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetConfigArguments)
include(GetGitVersion)

set(OMMP_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(OMMP_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
get_git_version(OMMP_VERSION OMMP_COMMIT)

project(
  "openmmpol"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "${OMMP_VERSION}"
  DESCRIPTION "Modern Fortran library for polarizable QM/MM simulations, commit ${OMMP_COMMIT}"
)

add_subdirectory("config")
get_config_arguments()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/openmmpol_const.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/openmmpol_const.h"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/version.h"
  @ONLY
)

file(COPY ${OMMP_ROOT_DIR}/include/openmmpol.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/)

if(NOT TARGET "OpenMP::OpenMP_Fortran")
  find_package("OpenMP" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

if(NOT TARGET "HDF5::HDF5")
  find_package("HDF5" REQUIRED)
endif()

if(NOT TARGET "cJSON")
  find_package("cJSON" REQUIRED)
endif()

if(NOT TARGET "OpenSSL")
  find_package("OpenSSL" REQUIRED)
endif()

# find_package("Python" REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "HDF5::HDF5"
  "cJSON"
  "OpenSSL"
)

# Collect sources of the project
set(srcs)
add_subdirectory("src")

# Setup main target, the openmmpol library
add_library(
  "${PROJECT_NAME}"
  "${srcs}"
)

set_target_properties(
  "${PROJECT_NAME}"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include"
)

# TODO: ugly, refactor
# if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018 -fno-realloc-lhs -fall-intrinsics")
#   target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic) # Very strict check on code
#   #target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2003 -fno-realloc-lhs) # Very strict check on code
#   target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-g -Wconversion>)
#   if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#       set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Warray-temporaries -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ffree-line-length-0 -fcheck=all")
#   endif()
#   # target_link_libraries(openmmpol stdc++)
#   set (lib-deps 
#       "${lib-deps}"
#       "stdc++")
# endif()

target_link_libraries(
  "${PROJECT_NAME}"
  PUBLIC
  "${lib-deps}"
)

target_include_directories(
  "${PROJECT_NAME}"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

target_include_directories(
  "${PROJECT_NAME}"
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src/
)

install(
  TARGETS
  "${PROJECT_NAME}"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)

# Python interface
if(WITH_PYTHON)
    message("-- Building with Python interface")
    add_subdirectory("utils/pythonapi")
endif()

# Tests
if(WITH_TESTS)
    message("-- Building with tests")
    enable_testing()
    add_subdirectory("tests")
endif()

