cmake_minimum_required(VERSION 3.18...3.29)

# Load modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(BuildUtils)
include(GetGitVersion)

# Get git version
set(OMMP_ROOT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
get_git_version(OMMP_VERSION OMMP_COMMIT)

# Define project
project(
  openmmpol
  LANGUAGES Fortran C CXX
  VERSION "${OMMP_VERSION}"
  DESCRIPTION "Modern Fortran library for polarizable QM/MM simulations, version commit: ${OMMP_COMMIT}"
)
include(GNUInstallDirs)

add_subdirectory("config")
get_config_arguments()
setup_build_type()
load_toolchain_settings()
setup_global_compiler_flags()


# /// Header file manipulation 
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
  "cJSON::cJSON"
  "OpenSSL::SSL"
)

if(NOT TARGET "OpenMP::OpenMP_Fortran")
  find_package("OpenMP" REQUIRED)
endif()

if(NOT TARGET "LAPACK::LAPACK")
  find_package("LAPACK" REQUIRED)
endif()

if(WITH_HDF5)
  list(APPEND lib-deps "hdf5::hdf5_fortran")

  if(NOT TARGET "hdf5::hdf5_fortran")
    find_package(HDF5 COMPONENTS Fortran)
  endif()

endif()

if(NOT TARGET "cJSON::cJSON")
  find_package("cJSON" REQUIRED)
endif()

if(NOT TARGET "OpenSSL::SSL")
  find_package("OpenSSL" REQUIRED)
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Collect sources of the project
set(srcs)
add_subdirectory("src")

# Preprocess all sources
# (mostly needed for Linux w\ makefile generator)
set_source_files_properties(
  ${srcs}
  PROPERTIES Fortran_PREPROCESS ON
)

# Setup main target, the openmmpol library
add_library(
  "${PROJECT_NAME}"
  SHARED
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

target_include_directories(
  "${PROJECT_NAME}"
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

target_link_libraries(
  "${PROJECT_NAME}"
  PUBLIC
  "${lib-deps}"
)

# Add executables
add_subdirectory("app")

install(
  TARGETS
  "${PROJECT_NAME}"
  EXPORT
  "${PROJECT_NAME}-targets"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/${module-dir}"
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}/${module-dir}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/${module-dir}"
)

install(
  EXPORT
  "${PROJECT_NAME}-targets"
  NAMESPACE
  "${PROJECT_NAME}::"
  # DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${module-dir}"
)

install(
  DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}"
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/app/ommp_pp"
  "${CMAKE_CURRENT_BINARY_DIR}/app/ommp_xyz2mmp"
  DESTINATION "${CMAKE_INSTALL_BINDIR}/${module-dir}"
  PERMISSIONS WORLD_EXECUTE
)

# Python interface
if(WITH_PYTHON)
    message("-- Building with Python interface")
    add_subdirectory("pythonapi")
endif()

# Tests
if(WITH_TESTS)
    message("-- Building with tests")
    enable_testing()
    # add_subdirectory("tests")
    include(tests/CMakeLists.txt)
endif()

