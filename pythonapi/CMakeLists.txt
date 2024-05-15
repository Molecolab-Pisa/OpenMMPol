find_package(pybind11 CONFIG REQUIRED)
pybind11_add_module("py${PROJECT_NAME}"
                     MODULE
                     "${CMAKE_CURRENT_SOURCE_DIR}/src/pyommp_interface.cpp" 
                     )
target_link_libraries("py${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}")

set_target_properties(
  "py${PROJECT_NAME}"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "py${PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
)

target_include_directories(
  "py${PROJECT_NAME}"
  PUBLIC
  $<BUILD_INTERFACE:${OMMP_ROOT_DIR}/include>
  $<BUILD_INTERFACE:${OMMP_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${module-dir}>
)

set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} -m pip install --no-deps --prefix ${prefix} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")