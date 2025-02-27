find_package(pybind11 CONFIG REQUIRED)
add_subdirectory("pyopenmmpol")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyopenmmpol/ommp_ip DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pyopenmmpol/)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/pyopenmmpol/topology_prm_assign DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/pyopenmmpol/)
set(prefix ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_PREFIX})
set(cmake-command "
  execute_process(
    COMMAND ${Python3_EXECUTABLE} -m pip install .
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
")
add_custom_target(PythonPackage
                  COMMAND ${Python3_EXECUTABLE} -m build ${CMAKE_CURRENT_BINARY_DIR}
                  DEPENDS "py${PROJECT_NAME}"
                  COMMENT "Generate Python Dist Files"
                 )
install(CODE "${cmake-command}"
        COMPONENT Python)
