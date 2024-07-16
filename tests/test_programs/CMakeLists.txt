# Add executable targets
add_executable(C_test_SI_init "tests/test_programs/C/test_SI_init.c")
add_executable(C_test_SI_potential "tests/test_programs/C/test_SI_potential.c")
add_executable(C_test_SI_geomgrad "tests/test_programs/C/test_SI_geomgrad.c")
add_executable(C_test_SI_geomgrad_num "tests/test_programs/C/test_SI_geomgrad_num.c")

# Link all executables to openmmpol
target_link_libraries(C_test_SI_init openmmpol)
target_link_libraries(C_test_SI_potential openmmpol)
target_link_libraries(C_test_SI_geomgrad openmmpol)
target_link_libraries(C_test_SI_geomgrad_num openmmpol)

# Put all targets into a proper directory
set_target_properties(C_test_SI_init
                    C_test_SI_potential
                    C_test_SI_geomgrad
                    C_test_SI_geomgrad_num
                    PROPERTIES
                    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# CI custom targets; do not touch
add_custom_target(C_test_programs DEPENDS C_test_SI_init
                                          C_test_SI_potential
                                          C_test_SI_geomgrad
                                          C_test_SI_geomgrad_num)


# Add executable targets
add_executable(F03_test_SI_init "tests/test_programs/F03/test_SI_init.f90")
add_executable(F03_test_SI_potential "tests/test_programs/F03/test_SI_potential.f90")
add_executable(F03_test_SI_geomgrad "tests/test_programs/F03/test_SI_geomgrad.f90")
add_executable(F03_test_SI_geomgrad_num "tests/test_programs/F03/test_SI_geomgrad_num.f90")

# Link all executables to openmmpol
target_link_libraries(F03_test_SI_init openmmpol)
target_link_libraries(F03_test_SI_potential openmmpol)
target_link_libraries(F03_test_SI_geomgrad openmmpol)
target_link_libraries(F03_test_SI_geomgrad_num openmmpol)

# Put all targets into a proper directory
set_target_properties(F03_test_SI_init
                      F03_test_SI_potential
                      F03_test_SI_geomgrad
                      F03_test_SI_geomgrad_num
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# CI custom targets; do not touch
add_custom_target(F03_test_programs DEPENDS F03_test_SI_init
                                            F03_test_SI_potential
                                            F03_test_SI_geomgrad
                                            F03_test_SI_geomgrad_num)
