## Download of FF parameters from Tinker repository
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoebabio18.prm
     amoebabio18.prm
     EXPECTED_HASH MD5=490ba93448bd71f14398d74f9cb8a523) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoeba09.prm
     amoeba09.prm
     EXPECTED_HASH MD5=91eaeec4edd2c581863c6cd4e0b2d975) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amber99sb.prm
     amber98sb.prm
     EXPECTED_HASH MD5=e9330e340a3f4045cc8025dc9ce42e33) 

#==== MEMORY TEST ====#
set(VALGRINDOPT --track-origins=yes --leak-check=full --show-leak-kinds=all) 
set(TESTLANG "C" CACHE STRING "Language to be used for the test suite, the very same test programs are implemented in C99 (C) and Fortran2003 (F03).")
add_test(NAME memcheck_NMA_loading_AMOEBA
         COMMAND valgrind ${VALGRINDOPT} bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp Testing/output)

add_test(NAME memcheck_NMA_loading_WANG
    COMMAND valgrind ${VALGRINDOPT} bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp Testing/output)

add_test(NAME memcheck_NMA_ene_full_AMOEBA
         COMMAND valgrind ${VALGRINDOPT} bin/${TESTLANG}_test_potential_xyz
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
         amoeba09.prm
         Testing/output)
add_test(NAME memcheck_NMA_ipd1_AMOEBA_CG
         COMMAND valgrind ${VALGRINDOPT} bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt)

#==== REGRESSION TEST ====#
# Tests on parameters loading
set(RTOL 1e-5)
set(ATOL 1e-6)

add_test(NMA_loading_AMOEBA
    bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp Testing/output)
add_test(NAME NMA_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_AMOEBA.ref)
add_test(NMA_loading_WANG_AL
    bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp Testing/output)
add_test(NAME NMA_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_WANG_AL.ref)

add_test(1crn_loading_AMOEBA 
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp Testing/output)
add_test(NAME 1crn_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_AMOEBA.ref)
add_test(1crn_loading_WANG_AL 
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp Testing/output)
add_test(NAME 1crn_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_WANG_AL.ref)

add_test(1ubq_loading_AMOEBA 
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp Testing/output)
add_test(NAME 1ubq_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_AMOEBA.ref)
add_test(1ubq_loading_WANG_AL 
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp Testing/output)
add_test(NAME 1ubq_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_WANG_AL.ref)

add_test(1ao6_loading_AMOEBA
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA.mmp Testing/output)
add_test(NAME 1ao6_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA.ref)
add_test(1ao6_loading_WANG_AL 
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL.mmp Testing/output)
add_test(NAME 1ao6_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_WANG_AL.ref)
add_test(1ao6_loading_AMOEBA_CUT10
         bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp Testing/output)
add_test(NAME 1ao6_loading_AMOEBA_CUT10_comp 
    COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA_CUT10.ref)
add_test(1ao6_loading_WANG_AL_CUT10
    bin/${TESTLANG}_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp Testing/output)
add_test(NAME 1ao6_loading_WANG_AL_CUT10_comp 
    COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_WANG_AL_CUT10.ref)

# Tests on loading from xyz/prm files

add_test(NMA_xyzloading_AMOEBA
    bin/${TESTLANG}_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz amoeba09.prm Testing/output)
set_tests_properties(NMA_xyzloading_AMOEBA PROPERTIES DEPENDS amoeba09.prm)
add_test(NAME NMA_xyzloading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_AMOEBA.ref)

add_test(1crn_xyzloading_AMOEBA
    bin/${TESTLANG}_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz amoebabio18.prm Testing/output)
set_tests_properties(1crn_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
add_test(NAME 1crn_xyzloading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_AMOEBA.ref)

add_test(1ubq_xyzloading_AMOEBA
    bin/${TESTLANG}_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz amoebabio18.prm Testing/output)
set_tests_properties(1ubq_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
add_test(NAME 1ubq_xyzloading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_AMOEBA.ref)

add_test(1ao6_xyzloading_AMOEBA
    bin/${TESTLANG}_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1ao6/input.xyz amoebabio18.prm Testing/output)
set_tests_properties(1ao6_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
add_test(NAME 1ao6_xyzloading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA.ref)

# Tests on calculation of IPD
add_test(NAME NMA_ipd0_AMOEBA_CG
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_AMOEBA_CG_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_CG
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_CG_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_INVERSION
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output inversion)
add_test(NAME NMA_ipd0_AMOEBA_INVERSION_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_INVERSION
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output inversion
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_INVERSION_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_DIIS
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output diis)
add_test(NAME NMA_ipd0_AMOEBA_DIIS_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_DIIS
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output diis
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_DIIS_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_CUT
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA_CUT.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_AMOEBA_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA_CUT.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_CUT
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA_CUT.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA_CUT.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt)
add_test(NAME NMA_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_WANG_AL_CUT
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL_CUT.mmp 
         Testing/output cg)
     add_test(NAME NMA_ipd0_WANG_AL_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_WANG_AL_CUT.ref 
         ${RTOL} ${ATOL})
     add_test(NAME NMA_ipd1_WANG_AL_CUT
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL_CUT.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
     add_test(NAME NMA_ipd1_WANG_AL_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_WANG_AL_CUT.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME 1crn_ipd0_AMOEBA
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME 1crn_ipd0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ipd1_AMOEBA
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ipd1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_1_AMOEBA.ref
         ${RTOL} ${ATOL})

add_test(NAME 1crn_ipd0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME 1crn_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ipd1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt)
add_test(NAME 1crn_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ubq_ipd0_AMOEBA
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME 1ubq_ipd0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ipd1_AMOEBA
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ipd1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_1_AMOEBA.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ubq_ipd0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME 1ubq_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ipd1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt)
add_test(NAME 1ubq_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ao6_ipd0_WANG_AL_CUT10
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output cg)
add_test(NAME 1ao6_ipd0_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_0_WANG_AL_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd1_WANG_AL_CUT10
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt)
add_test(NAME 1ao6_ipd1_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_1_WANG_AL_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd0_AMOEBA_CUT10
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output cg)
add_test(NAME 1ao6_ipd0_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_0_AMOEBA_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd1_AMOEBA_CUT10
         COMMAND  bin/${TESTLANG}_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt)
add_test(NAME 1ao6_ipd1_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_1_AMOEBA_CUT10.ref
         ${RTOL} ${ATOL})

# Tests on calculation of energy
add_test(NAME NMA_ene0_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME NMA_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene1_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME NMA_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene0_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME 1crn_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene1_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME 1crn_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene0_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME 1ubq_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene1_AMOEBA
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene0_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME 1ubq_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene1_WANG_AL
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME 1ao6_ene0_AMOEBA_CUT10
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output)
add_test(NAME 1ao6_ene0_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_0_AMOEBA_CUT10.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene1_AMOEBA_CUT10
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt) 
add_test(NAME 1ao6_ene1_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_1_AMOEBA_CUT10.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene0_WANG_AL_CUT10
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output)
add_test(NAME 1ao6_ene0_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_0_WANG_AL_CUT10.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene1_WANG_AL_CUT10
         COMMAND  bin/${TESTLANG}_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt) 
add_test(NAME 1ao6_ene1_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_1_WANG_AL_CUT10.ref
         Testing/output 
         ${RTOL} ${ATOL})

# Tests on the full potential
 add_test(NAME NMA_ene_full_AMOEBA
          COMMAND  bin/${TESTLANG}_test_potential_xyz
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
          amoeba09.prm
          Testing/output)
 add_test(NAME NMA_ene_full_AMOEBA_comp_eb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EB
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ev
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EV
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_em
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EM
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ep
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EP
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_eub
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EUB
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ea
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EA
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_eba
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EBA
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_eopb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EOPB
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ept
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EPT
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_et
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          ET
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ett
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          ETT
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_eat
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EAT
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_ebt
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          EBT
          ${RTOL} ${ATOL})
 add_test(NAME NMA_ene_full_AMOEBA_comp_etot
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
          Testing/output 
          ETOT
          ${RTOL} ${ATOL})
 
 add_test(NAME 1crn_ene_full_AMOEBA
          COMMAND  bin/${TESTLANG}_test_potential_xyz
          ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz
          amoebabio18.prm
          Testing/output)
 add_test(NAME 1crn_ene_full_AMOEBA_comp_eb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EB
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ev
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EV
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_em
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EM
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ep
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EP
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_eub
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EUB
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ea
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EA
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_eba
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EBA
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_eopb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EOPB
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ept
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EPT
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_et
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          ET
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ett
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          ETT
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_eat
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EAT
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_ebt
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          EBT
          ${RTOL} ${ATOL})
 add_test(NAME 1crn_ene_full_AMOEBA_comp_etot
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
          Testing/output 
          ETOT
          2.0e-5 3.0e-4)
 
 add_test(NAME 1ubq_ene_full_AMOEBA
          COMMAND  bin/${TESTLANG}_test_potential_xyz
          ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz
          amoebabio18.prm
          Testing/output)
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_eb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EB
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ev
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EV
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_em
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EM
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ep
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EP
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_eub
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EUB
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ea
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EA
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_eba
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EBA
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_eopb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EOPB
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ept
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EPT
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_et
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          ET
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ett
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          ETT
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_eat
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EAT
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_ebt
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          EBT
          ${RTOL} ${ATOL})
 add_test(NAME 1ubq_ene_full_AMOEBA_comp_etot
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
          Testing/output 
          ETOT
          ${RTOL} ${ATOL})
 
 add_test(NAME 2d94_ene_full_AMOEBA
          COMMAND  bin/${TESTLANG}_test_potential_xyz
          ${CMAKE_SOURCE_DIR}/tests/2d94/input.xyz
          amoebabio18.prm
          Testing/output)
 add_test(NAME 2d94_ene_full_AMOEBA_comp_eb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EB
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ev
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EV
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_em
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EM
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ep
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EP
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_eub
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EUB
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ea
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EA
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_eba
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EBA
          1e-4 1e-5)
 add_test(NAME 2d94_ene_full_AMOEBA_comp_eopb
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EOPB
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ept
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EPT
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_et
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          ET
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ett
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          ETT
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_eat
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EAT
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_ebt
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          EBT
          ${RTOL} ${ATOL})
 add_test(NAME 2d94_ene_full_AMOEBA_comp_etot
          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
          Testing/output 
          ETOT
          ${RTOL} ${ATOL})
 
 # Tests on the full potential (after reloading from HDF5)
 if (HDF5_WORKS)
     add_test(NAME NMA_ene_full_AMOEBA_hdf5_dump
              COMMAND  bin/ommp_pp
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
              amoeba09.prm
              Testing/NMA.hdf5)
     add_test(NAME NMA_ene_full_AMOEBA_hdf5_load
              COMMAND  bin/${TESTLANG}_test_potential_hdf5
              Testing/NMA.hdf5
              Testing/output_hdf5)
     add_test(NAME NMA_ene_full_AMOEBA_comp_eb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5
              EB
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ev_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EV
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_em_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EM
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ep_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EP
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_eub_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EUB
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ea_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EA
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_eba_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBA
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_eopb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EOPB
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ept_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EPT
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_et_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ET
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ett_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETT
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_eat_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EAT
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_ebt_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBT
              ${RTOL} ${ATOL})
     add_test(NAME NMA_ene_full_AMOEBA_comp_etot_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETOT
              ${RTOL} ${ATOL})
 
     add_test(NAME 1crn_ene_full_AMOEBA_hdf5_dump
              COMMAND  bin/ommp_pp
              ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz
              amoebabio18.prm
              Testing/1crn.hdf5)
     add_test(NAME 1crn_ene_full_AMOEBA_hdf5_load
              COMMAND  bin/${TESTLANG}_test_potential_hdf5
              Testing/1crn.hdf5
              Testing/output_hdf5)
     add_test(NAME 1crn_ene_full_AMOEBA_comp_eb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EB
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ev_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EV
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_em_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EM
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ep_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EP
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_eub_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EUB
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ea_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EA
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_eba_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBA
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_eopb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EOPB
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ept_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EPT
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_et_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ET
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ett_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETT
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_eat_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EAT
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_ebt_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBT
              ${RTOL} ${ATOL})
     add_test(NAME 1crn_ene_full_AMOEBA_comp_etot_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETOT
              2.0e-5 3.0e-4)
 
     add_test(NAME 1ubq_ene_full_AMOEBA_hdf5_dump
              COMMAND  bin/ommp_pp
              ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz
              amoebabio18.prm
              Testing/1ubq.hdf5)
     add_test(NAME 1ubq_ene_full_AMOEBA_hdf5_load
              COMMAND  bin/${TESTLANG}_test_potential_hdf5
              Testing/1ubq.hdf5
              Testing/output_hdf5)
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EB
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ev_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EV
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_em_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EM
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ep_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EP
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eub_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EUB
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ea_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EA
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eba_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBA
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eopb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EOPB
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ept_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EPT
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_et_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ET
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ett_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETT
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eat_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EAT
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ebt_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBT
              ${RTOL} ${ATOL})
     add_test(NAME 1ubq_ene_full_AMOEBA_comp_etot_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETOT
              ${RTOL} ${ATOL})
 
     add_test(NAME 2d94_ene_full_AMOEBA_hdf5_dump
              COMMAND  bin/ommp_pp
              ${CMAKE_SOURCE_DIR}/tests/2d94/input.xyz
              amoebabio18.prm
              Testing/2d94.hdf5)
     add_test(NAME 2d94_ene_full_AMOEBA_hdf5_load
              COMMAND  bin/${TESTLANG}_test_potential_hdf5
              Testing/2d94.hdf5
              Testing/output_hdf5)
     add_test(NAME 2d94_ene_full_AMOEBA_comp_eb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EB
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ev_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EV
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_em_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EM
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ep_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EP
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_eub_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EUB
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ea_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EA
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_eba_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBA
              1e-4 1e-5)
     add_test(NAME 2d94_ene_full_AMOEBA_comp_eopb_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EOPB
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ept_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EPT
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_et_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ET
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ett_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETT
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_eat_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EAT
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_ebt_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              EBT
              ${RTOL} ${ATOL})
     add_test(NAME 2d94_ene_full_AMOEBA_comp_etot_hdf5
              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
              Testing/output_hdf5 
              ETOT
              ${RTOL} ${ATOL})
endif ()


# Test on Gradients
set(RTOL 1e-3)
set(ATOL 1e-4)
add_test(NAME NMA_geomgradn_WANG_AL
         COMMAND  bin/${TESTLANG}_test_geomgrad_num
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME NMA_geomgradn_WANG_AL_CUT
         COMMAND  bin/${TESTLANG}_test_geomgrad_num
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL_CUT.mmp 
         Testing/output)
add_test(NAME NMA_geomgradn_AMOEBA
         COMMAND  bin/${TESTLANG}_test_geomgrad_num
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME NMA_geomgradn_AMOEBA_CUT
         COMMAND  bin/${TESTLANG}_test_geomgrad_num
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA_CUT.mmp 
         Testing/output)
add_test(NAME NMA_geomgradn_full_AMOEBA
          COMMAND  bin/${TESTLANG}_test_geomgrad_xyz_num
          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
          amoeba09.prm
          Testing/output)
 
add_test(NAME NMA_geomgrad_full_AMOEBA
         COMMAND  bin/${TESTLANG}_test_geomgrad_xyz
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
         amoeba09.prm
         Testing/output)
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_eb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EB
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ev
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EV
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_em
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EM
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ep
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EP
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_eub
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EUB
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ea
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EA
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_eba
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EBA
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_eopb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EOPB
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ept
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EPT
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_et
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         ET
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ett
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         ETT
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_eat
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EAT
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_ebt
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         EBT
         ${RTOL} ${ATOL})
add_test(NAME NMA_geomgrad_full_AMOEBA_comp_etot
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
         Testing/output 
         ETOT
         ${RTOL} ${ATOL})
 
add_test(NAME 1crn_geomgrad_full_AMOEBA
         COMMAND  bin/${TESTLANG}_test_geomgrad_xyz
         ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz
         amoebabio18.prm
         Testing/output)
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_eb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EB
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ev
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EV
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_em
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EM
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ep
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EP
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_eub
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EUB
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ea
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EA
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_eba
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EBA
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_eopb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EOPB
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ept
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EPT
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_et
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         ET
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ett
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         ETT
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_eat
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EAT
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_ebt
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         EBT
         ${RTOL} ${ATOL})
add_test(NAME 1crn_geomgrad_full_AMOEBA_comp_etot
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
         Testing/output 
         ETOT
         ${RTOL} ${ATOL})
 
add_test(NAME 1ubq_geomgrad_full_AMOEBA
         COMMAND  bin/${TESTLANG}_test_geomgrad_xyz
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz
         amoebabio18.prm
         Testing/output)
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_eb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EB
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ev
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EV
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_em
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EM
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ep
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EP
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_eub
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EUB
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ea
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EA
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_eba
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EBA
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_eopb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EOPB
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ept
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EPT
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_et
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         ET
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ett
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         ETT
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_eat
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EAT
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_ebt
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         EBT
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_geomgrad_full_AMOEBA_comp_etot
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
         Testing/output 
         ETOT
         ${RTOL} ${ATOL})
 
add_test(NAME 2d94_geomgrad_full_AMOEBA
         COMMAND  bin/${TESTLANG}_test_geomgrad_xyz
         ${CMAKE_SOURCE_DIR}/tests/2d94/input.xyz
         amoebabio18.prm
         Testing/output)
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_eb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EB
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ev
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EV
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_em
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EM
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ep
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EP
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_eub
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EUB
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ea
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EA
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_eba
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EBA
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_eopb
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EOPB
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ept
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EPT
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_et
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         ET
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ett
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         ETT
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_eat
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EAT
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_ebt
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         EBT
         ${RTOL} ${ATOL})
add_test(NAME 2d94_geomgrad_full_AMOEBA_comp_etot
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_geomgrad.py
         ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
         Testing/output 
         ETOT
         ${RTOL} ${ATOL})
