find_package(pybind11 CONFIG REQUIRED)
add_library("py${PROJECT_NAME}" PUBLIC)
target_link_libraries("py${PROJECT_NAME}" INTERFACE "${PROJECT_NAME}")
pybind11_add_module("py${PROJECT_NAME}" MODULE "${CMAKE_CURRENT_SOURCE_DIR}/src/pyommp_interface.cpp")

set_target_properties(
  "py${PROJECT_NAME}"
  PROPERTIES
  POSITION_INDEPENDENT_CODE TRUE
  OUTPUT_NAME "$py{PROJECT_NAME}"
  VERSION "${PROJECT_VERSION}"
  SOVERSION "${PROJECT_VERSION_MAJOR}"
)

set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${PYTHON_INTERPRETER} -m pip install --no-deps --prefix ${prefix} .
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
")

install(CODE "${cmake-command}")