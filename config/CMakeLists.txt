if(NOT DEFINED "${PROJECT_NAME}-dependency-method")
  set(
    "${PROJECT_NAME}-dependency-method"
    "subproject" "cmake" "pkgconf" "fetch"
  )
endif()

set(
  "${PROJECT_NAME}-module-dir"
  "${PROJECT_NAME}/${PROJECT_NAME}-${OMMP_VERSION}_${CMAKE_Fortran_COMPILER_ID}-${CMAKE_Fortran_COMPILER_VERSION}"
  CACHE STRING
  "Subdirectory to install generated module files to"
)
set(
  module-dir
  "${${PROJECT_NAME}-module-dir}"
)
set(module-dir "${module-dir}" PARENT_SCOPE)

# Set build type as CMake does not provide defaults
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(
    CMAKE_BUILD_TYPE "RelWithDebInfo"
    CACHE STRING "Build type to be used."
    FORCE
  )
  message(
    STATUS
    "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified."
  )
endif()

include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.pc"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.pc"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY SameMinorVersion
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/template.pc"
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  @ONLY
)

install(
  FILES
  "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
)

