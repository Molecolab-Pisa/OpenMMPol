if(${TESTLANG} STREQUAL "C")

add_executable(${TESTLANG}_test_SI_init "tests/test_programs/C/test_SI_init.c")
add_executable(${TESTLANG}_test_SI_potential "tests/test_programs/C/test_SI_potential.c")
add_executable(${TESTLANG}_test_SI_geomgrad "tests/test_programs/C/test_SI_geomgrad.c")
add_executable(${TESTLANG}_test_SI_geomgrad_num "tests/test_programs/C/test_SI_geomgrad_num.c")

elseif(${TESTLANG} STREQUAL "F03")

add_executable(${TESTLANG}_test_SI_init "tests/test_programs/F03/test_SI_init.f90")
add_executable(${TESTLANG}_test_SI_potential "tests/test_programs/F03/test_SI_potential.f90")
add_executable(${TESTLANG}_test_SI_geomgrad "tests/test_programs/F03/test_SI_geomgrad.f90")
add_executable(${TESTLANG}_test_SI_geomgrad_num "tests/test_programs/F03/test_SI_geomgrad_num.f90")

endif()

target_link_libraries(${TESTLANG}_test_SI_init openmmpol)
target_link_libraries(${TESTLANG}_test_SI_potential openmmpol)
target_link_libraries(${TESTLANG}_test_SI_geomgrad openmmpol)
target_link_libraries(${TESTLANG}_test_SI_geomgrad_num openmmpol)

set_target_properties(${TESTLANG}_test_SI_init
                      ${TESTLANG}_test_SI_potential
                      ${TESTLANG}_test_SI_geomgrad
                      ${TESTLANG}_test_SI_geomgrad_num
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)