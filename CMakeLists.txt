# The whole cmake structure of the project is
# heavily borrowing from
# https://github.com/dftd4/dftd4/

cmake_minimum_required(VERSION 3.18)

# Get package version from Git

# TODO: get_git_version does not work, fix it
set(OMMP_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
include("${CMAKE_CURRENT_SOURCE_DIR}/config/cmake/GetGitVersion.cmake")
get_git_version(OMMP_VERSION)

project(
  "openmmpol"
  LANGUAGES "Fortran" "C" "CXX"
  VERSION "${OMMP_VERSION}"
  DESCRIPTION "Modern Fortran library for polarizable QM/MM simulations"
)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/config/")

# Get config arguments
get_config_arguments()

if(NOT TARGET "OpenMP::OpenMP_Fortran")
  find_package("OpenMP" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

if(NOT TARGET "HDF5::HDF5")
  find_package("HDF5" REQUIRED)
endif()

if(NOT TARGET "cJSON::cJSON")
  find_package("cJSON" REQUIRED)
endif()

if(NOT TARGET "OpenSSL::OpenSSL")
  find_package("OpenSSL" REQUIRED)
endif()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "HDF5::HDF5"
  "cJSON::cJSON"
  "OpenSSL::OpenSSL"
)

# Python interpeter
set(PYTHON_INTERPRETER "python3" CACHE STRING
  "Python interpreter to use for installing and testing python components")

# Collect sources of the project
set(srcs)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/")

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
    add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/utils/pythonapi")
endif()

# Tests
if(WITH_TESTS)
    enable_testing()
    add_subdirectory("tests")
endif()

