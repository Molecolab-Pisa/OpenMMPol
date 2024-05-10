@PACKAGE_INIT@

set(
  "@PROJECT_NAME@_INCLUDE_DIRS"
  "@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_INCLUDEDIR@"
  "@CMAKE_INSTALL_PREFIX@/@CMAKE_INSTALL_INCLUDEDIR@/@module-dir@"
)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")

  include(CMakeFindDependencyMacro)

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
  
endif()