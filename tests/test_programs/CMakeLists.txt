foreach(lang "C" "F03")
    if(${lang} STREQUAL "C")
        set(ext "c")
    elseif(${lang} STREQUAL "F03")
        set(ext "f90")
    endif()

    add_executable(${lang}_test_SI_init "tests/test_programs/${lang}/test_SI_init.${ext}")
    add_executable(${lang}_test_SI_potential "tests/test_programs/${lang}/test_SI_potential.${ext}")
    add_executable(${lang}_test_SI_geomgrad "tests/test_programs/${lang}/test_SI_geomgrad.${ext}")
    add_executable(${lang}_test_SI_geomgrad_num "tests/test_programs/${lang}/test_SI_geomgrad_num.${ext}")
    
    target_link_libraries(${lang}_test_SI_init openmmpol)
    target_link_libraries(${lang}_test_SI_potential openmmpol)
    target_link_libraries(${lang}_test_SI_geomgrad openmmpol)
    target_link_libraries(${lang}_test_SI_geomgrad_num openmmpol)

    set_target_properties(${lang}_test_SI_init
                        ${lang}_test_SI_potential
                        ${lang}_test_SI_geomgrad
                        ${lang}_test_SI_geomgrad_num
                        PROPERTIES
                        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

    add_custom_target(${lang}_test_programs DEPENDS ${lang}_test_SI_init
                                                        ${lang}_test_SI_potential
                                                        ${lang}_test_SI_geomgrad
                                                        ${lang}_test_SI_geomgrad_num)
endforeach()
