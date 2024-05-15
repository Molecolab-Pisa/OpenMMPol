# The whole cmake structure of this project is
# heavily borrowing from
# https://github.com/dftd4/dftd4/
# and 
# https://github.com/dftbplus/dftbplus/

cmake_minimum_required(VERSION 3.18...3.29)
# cmake_minimum_required(VERSION 3.27)

# if(${CMAKE_VERSION} VERSION_LESS 3.12)
    # cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
# endif()

# Load modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(BuildUtils)
include(GetGitVersion)

set(OMMP_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
# set(OMMP_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")
get_git_version(OMMP_VERSION OMMP_COMMIT)

project(
  openmmpol
  LANGUAGES Fortran C CXX
  VERSION "${OMMP_VERSION}"
  DESCRIPTION "Modern Fortran library for polarizable QM/MM simulations, version commit: ${OMMP_COMMIT}"
)

add_subdirectory("config")
get_config_arguments()
setup_build_type()
load_toolchain_settings()
setup_global_compiler_flags()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/openmmpol_const.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/openmmpol_const.h"
  @ONLY
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/version.h"
  @ONLY
)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/openmmpol.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/include/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/include/f_cart_components.h DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "cJSON"
  "OpenSSL"
  # TODO: add generator expression for HDF5
  # that depends on WITH_HDF5 variable?
)

if(NOT TARGET "OpenMP::OpenMP_Fortran")
  find_package("OpenMP" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

if(WITH_HDF5)
  list(APPEND lib-deps "HDF5::HDF5")

  if(NOT TARGET "HDF5::HDF5")
    find_package("HDF5" REQUIRED)
  endif()

endif()

if(NOT TARGET "cJSON")
  find_package("cJSON" REQUIRED)
endif()

if(NOT TARGET "OpenSSL")
  find_package("OpenSSL" REQUIRED)
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Collect sources of the project
set(srcs)
# set(srcs_preprocess)
add_subdirectory("src")

set_source_files_properties(
  # ${srcs_preprocess}
  ${srcs}
  PROPERTIES Fortran_PREPROCESS ON
)

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
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod"
)

target_link_libraries(
  "${PROJECT_NAME}"
  PUBLIC
  "${lib-deps}"
)

target_include_directories(
  "${PROJECT_NAME}"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
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

