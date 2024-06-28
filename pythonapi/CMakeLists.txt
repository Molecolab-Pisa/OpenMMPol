find_package(pybind11 CONFIG REQUIRED)
add_subdirectory("pyopenmmpol")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${Python_EXECUTABLE} -m pip install .
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
")

install(CODE "${cmake-command}")