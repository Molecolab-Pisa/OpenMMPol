pybind11_add_module("py${PROJECT_NAME}"
                     MODULE
                     pyommp_interface.cpp
                     )

set_target_properties(
   "py${PROJECT_NAME}"
   PROPERTIES
   POSITION_INDEPENDENT_CODE TRUE
   OUTPUT_NAME "__py${PROJECT_NAME}"
   VERSION "${PROJECT_VERSION}"
   SOVERSION "${PROJECT_VERSION_MAJOR}"
 )
 
target_link_libraries("py${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}")

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
configure_file(__version__.py.in ${CMAKE_CURRENT_BINARY_DIR}/__version__.py @ONLY)