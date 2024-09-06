cmake_minimum_required(VERSION 3.18...3.29)
# Load modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(CMakePrintHelpers)
include(BuildUtils)
include(GetGitVersion)

# Get git version
get_git_version(OMMP_VERSION OMMP_COMMIT OMMP_VERSION_INTERNAL)
  
# Define project
project(
  openmmpol
  LANGUAGES Fortran C CXX
  VERSION ${OMMP_VERSION}
  DESCRIPTION "Modern Fortran library for polarizable QM/MM simulations, version commit: ${OMMP_COMMIT}"
)
include(GNUInstallDirs)

add_subdirectory("config")
get_config_arguments()

set(
  lib-deps
  "OpenMP::OpenMP_Fortran"
  "LAPACK::LAPACK"
  "cjson"
  "OpenSSL::SSL"
)

if(NOT TARGET OpenMP::OpenMP_Fortran)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

if(NOT TARGET LAPACK::LAPACK)
  find_package(LAPACK REQUIRED)
endif()

if(WITH_HDF5)
  list(APPEND lib-deps hdf5::hdf5_fortran)

  if(NOT TARGET hdf5::hdf5_fortran)
    find_package(HDF5 COMPONENTS Fortran)
  endif()

endif()

if(NOT TARGET cjson)
  find_package(cJSON REQUIRED)
endif()

if(NOT TARGET OpenSSL::SSL)
  find_package(OpenSSL REQUIRED)
endif()

find_package(Python COMPONENTS Interpreter Development REQUIRED)

setup_build_type()
load_toolchain_settings()
setup_global_compiler_flags()

# Collect sources of the project
set(srcs)
add_subdirectory("src")

# Explicitly require CMake to preprocess all sources;
# this is necessary for some generators
set_source_files_properties(
  ${srcs}
  PROPERTIES Fortran_PREPROCESS ON
)

# Setup main library target
add_library(
  ${PROJECT_NAME}
  SHARED
  "${srcs}"
)

set_target_properties(
  ${PROJECT_NAME}
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME ${PROJECT_NAME}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod
)

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/lib>
)

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
  "${lib-deps}"
)

# Add executables
add_subdirectory("app")

# Install targets 1
install(
  TARGETS
  ${PROJECT_NAME}
  EXPORT
  ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  INCLUDES DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# Install targets 2
install(
  EXPORT
  ${PROJECT_NAME}-targets
  NAMESPACE
  ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}
  COMPONENT Library
)

# Install CMake config files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/openmmpol-config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/openmmpol-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/${PROJECT_NAME}
        COMPONENT Library
        )

# Install the Fortran .mod
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/mod/ommp_interface.mod
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/
        COMPONENT Library)

# Install the headers
install(
  DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
  COMPONENT Library
)

# Install binaries
install(
  FILES
  ${CMAKE_CURRENT_BINARY_DIR}/app/ommp_xyz2mmp
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  PERMISSIONS WORLD_EXECUTE
  COMPONENT UtilityApp
)
if(WITH_HDF5)
    install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/app/ommp_pp
      DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
      PERMISSIONS WORLD_EXECUTE
      COMPONENT UtilityApp
    )
endif()

# Python interface
if(WITH_PYTHON)
    message("-- Building with Python interface")
    add_subdirectory("pythonapi")
endif()

# Tests
if(WITH_TESTS)
    message("-- Building with tests")
    enable_testing()
    # Currently, the test language is set here.
    # TODO: add an option to config.cmake to alternate between C and F03
    set(TESTLANG "C"
        CACHE STRING
        "Language to be used for the test suite,
        the very same test programs are implemented in C99 (C)
        and Fortran2003 (F03).")
    include(tests/CMakeLists.txt)
endif()

