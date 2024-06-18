include(CMakePackageConfigHelpers)

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

configure_file(openmmpol-config.cmake.in ${CMAKE_BINARY_DIR}/openmmpol-config.cmake @ONLY)

configure_file(
  openmmpol_const.h.in
  ${CMAKE_BINARY_DIR}/openmmpol_const.h
  @ONLY
)

configure_file(
  version.h.in
  ${CMAKE_BINARY_DIR}/version.h
  @ONLY
)

# Header file manipulation 
file(COPY ${CMAKE_SOURCE_DIR}/include/openmmpol.h DESTINATION ${CMAKE_BINARY_DIR}/include)
file(COPY ${CMAKE_SOURCE_DIR}/include/f_cart_components.h DESTINATION ${CMAKE_BINARY_DIR}/include)

write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake"
  VERSION "${PROJECT_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

