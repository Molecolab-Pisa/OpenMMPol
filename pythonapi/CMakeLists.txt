find_package(pybind11 CONFIG REQUIRED)
add_subdirectory("pyopenmmpol")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${Python_EXECUTABLE} -m pip install .
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
")
add_custom_target(PythonPackage
                  COMMAND ${Python_EXECUTABLE} -m build ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS "py${PROJECT_NAME}"
                  COMMENT "Generate Python Dist Files"
                 )
install(CODE "${cmake-command}"
        COMPONENT Python)
