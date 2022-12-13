## Download of FF parameters from Tinker repository
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoebabio18.prm
     amoebabio18.prm
     EXPECTED_HASH MD5=490ba93448bd71f14398d74f9cb8a523) 
file(DOWNLOAD https://raw.githubusercontent.com/TinkerTools/tinker/release/params/amoeba09.prm
     amoeba09.prm
     EXPECTED_HASH MD5=91eaeec4edd2c581863c6cd4e0b2d975) 

#==== MEMORY TEST ====#
set(VALGRINDOPT --track-origins=yes --leak-check=full --show-leak-kinds=all) 
add_test(NAME memcheck_NMA_loading_AMOEBA
         COMMAND valgrind ${VALGRINDOPT} bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp Testing/output)

add_test(NAME memcheck_NMA_loading_WANG
    COMMAND valgrind ${VALGRINDOPT} bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp Testing/output)

#TODO# add_test(NAME memcheck_NMA_ene_full_AMOEBA
#TODO#          COMMAND valgrind ${VALGRINDOPT} bin/C_test_potential_xyz
#TODO#          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
#TODO#          amoeba09.prm
#TODO#          Testing/output)
add_test(NAME memcheck_NMA_ipd1_AMOEBA_CG
         COMMAND valgrind ${VALGRINDOPT} bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt)

#==== REGRESSION TEST ====#
# Tests on parameters loading
set(RTOL 1e-5)
set(ATOL 1e-6)

add_test(NMA_loading_AMOEBA
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp Testing/output)
add_test(NAME NMA_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_AMOEBA.ref)
add_test(NMA_loading_WANG_AL
    bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp Testing/output)
add_test(NAME NMA_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_WANG_AL.ref)

add_test(1crn_loading_AMOEBA 
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp Testing/output)
add_test(NAME 1crn_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_AMOEBA.ref)
add_test(1crn_loading_WANG_AL 
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp Testing/output)
add_test(NAME 1crn_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_WANG_AL.ref)

add_test(1ubq_loading_AMOEBA 
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp Testing/output)
add_test(NAME 1ubq_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_AMOEBA.ref)
add_test(1ubq_loading_WANG_AL 
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp Testing/output)
add_test(NAME 1ubq_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_WANG_AL.ref)

add_test(1ao6_loading_AMOEBA
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA.mmp Testing/output)
add_test(NAME 1ao6_loading_AMOEBA_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA.ref)
add_test(1ao6_loading_WANG_AL 
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL.mmp Testing/output)
add_test(NAME 1ao6_loading_WANG_AL_comp 
         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_WANG_AL.ref)
add_test(1ao6_loading_AMOEBA_CUT10
         bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp Testing/output)
add_test(NAME 1ao6_loading_AMOEBA_CUT10_comp 
    COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA_CUT10.ref)
add_test(1ao6_loading_WANG_AL_CUT10
    bin/C_test_init ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp Testing/output)
add_test(NAME 1ao6_loading_WANG_AL_CUT10_comp 
    COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_WANG_AL_CUT10.ref)

# Tests on loading from xyz/prm files
#
#add_test(NMA_xyzloading_AMOEBA
#    bin/C_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz amoeba09.prm Testing/output)
#set_tests_properties(NMA_xyzloading_AMOEBA PROPERTIES DEPENDS amoeba09.prm)
#add_test(NAME NMA_xyzloading_AMOEBA_comp 
#         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/summary_AMOEBA.ref)
#
#add_test(1crn_xyzloading_AMOEBA
#    bin/C_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz amoebabio18.prm Testing/output)
#set_tests_properties(1crn_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
#add_test(NAME 1crn_xyzloading_AMOEBA_comp 
#         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1crn/summary_AMOEBA.ref)
#
#add_test(1ubq_xyzloading_AMOEBA
#    bin/C_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz amoebabio18.prm Testing/output)
#set_tests_properties(1ubq_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
#add_test(NAME 1ubq_xyzloading_AMOEBA_comp 
#         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ubq/summary_AMOEBA.ref)
#
#add_test(1ao6_xyzloading_AMOEBA
#    bin/C_test_init_xyz ${CMAKE_SOURCE_DIR}/tests/1ao6/input.xyz amoebabio18.prm Testing/output)
#set_tests_properties(1ao6_xyzloading_AMOEBA PROPERTIES DEPENDS amoebabio18.prm)
#add_test(NAME 1ao6_xyzloading_AMOEBA_comp 
#         COMMAND ${CMAKE_COMMAND} -E compare_files Testing/output ${CMAKE_SOURCE_DIR}/tests/1ao6/summary_AMOEBA.ref)
#
# Tests on calculation of IPD
add_test(NAME NMA_ipd0_AMOEBA_CG
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_AMOEBA_CG_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_CG
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_CG_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_INVERSION
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output inversion)
add_test(NAME NMA_ipd0_AMOEBA_INVERSION_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_INVERSION
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output inversion
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_INVERSION_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_DIIS
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output diis)
add_test(NAME NMA_ipd0_AMOEBA_DIIS_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_DIIS
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output diis
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_DIIS_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_AMOEBA_CUT
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA_CUT.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_AMOEBA_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_AMOEBA_CUT.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_AMOEBA_CUT
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA_CUT.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ipd1_AMOEBA_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_AMOEBA_CUT.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME NMA_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME NMA_ipd1_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt)
add_test(NAME NMA_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME NMA_ipd0_WANG_AL_CUT
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL_CUT.mmp 
         Testing/output cg)
     add_test(NAME NMA_ipd0_WANG_AL_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_0_WANG_AL_CUT.ref 
         ${RTOL} ${ATOL})
     add_test(NAME NMA_ipd1_WANG_AL_CUT
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL_CUT.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
     add_test(NAME NMA_ipd1_WANG_AL_CUT_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/IPD_1_WANG_AL_CUT.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME 1crn_ipd0_AMOEBA
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME 1crn_ipd0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ipd1_AMOEBA
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ipd1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_1_AMOEBA.ref
         ${RTOL} ${ATOL})

add_test(NAME 1crn_ipd0_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME 1crn_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ipd1_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt)
add_test(NAME 1crn_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ubq_ipd0_AMOEBA
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output cg)
add_test(NAME 1ubq_ipd0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ipd1_AMOEBA
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ipd1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_1_AMOEBA.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ubq_ipd0_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output cg)
add_test(NAME 1ubq_ipd0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_0_WANG_AL.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ipd1_WANG_AL
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt)
add_test(NAME 1ubq_ipd1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/IPD_1_WANG_AL.ref
         ${RTOL} ${ATOL})

add_test(NAME 1ao6_ipd0_WANG_AL_CUT10
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output cg)
add_test(NAME 1a06_ipd0_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_0_WANG_AL_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd1_WANG_AL_CUT10
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output cg
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt)
add_test(NAME 1ao6_ipd1_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_1_WANG_AL_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd0_AMOEBA_CUT10
         COMMAND  bin/C_test_ipd 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output cg)
add_test(NAME 1a06_ipd0_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/IPD_0_AMOEBA_CUT10.ref
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ipd1_AMOEBA_CUT10
         COMMAND  bin/C_test_ipd 
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
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME NMA_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene1_AMOEBA
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene0_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME NMA_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME NMA_ene1_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/EF_1.txt) 
add_test(NAME NMA_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene0_AMOEBA
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME 1crn_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene1_AMOEBA
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene0_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME 1crn_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1crn_ene1_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1crn/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1crn/EF_1.txt) 
add_test(NAME 1crn_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1crn/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene0_AMOEBA
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output)
add_test(NAME 1ubq_ene0_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_0_AMOEBA.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene1_AMOEBA
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_AMOEBA.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ene1_AMOEBA_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_1_AMOEBA.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene0_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output)
add_test(NAME 1ubq_ene0_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_0_WANG_AL.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ubq_ene1_WANG_AL
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ubq/input_WANG_AL.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ubq/EF_1.txt) 
add_test(NAME 1ubq_ene1_WANG_AL_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ubq/ENE_1_WANG_AL.ref
         Testing/output 
         ${RTOL} ${ATOL})

add_test(NAME 1ao6_ene0_AMOEBA_CUT10
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output)
add_test(NAME 1ao6_ene0_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_0_AMOEBA_CUT10.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene1_AMOEBA_CUT10
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_AMOEBA_CUT10.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt) 
add_test(NAME 1ao6_ene1_AMOEBA_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_1_AMOEBA_CUT10.ref
         Testing/output 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene0_WANG_AL_CUT10
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output)
add_test(NAME 1ao6_ene0_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         Testing/output 
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_0_WANG_AL_CUT10.ref 
         ${RTOL} ${ATOL})
add_test(NAME 1ao6_ene1_WANG_AL_CUT10
         COMMAND  bin/C_test_energy
         ${CMAKE_SOURCE_DIR}/tests/1ao6/input_WANG_AL_CUT10.mmp 
         Testing/output
         ${CMAKE_SOURCE_DIR}/tests/1ao6/EF_1.txt) 
add_test(NAME 1ao6_ene1_WANG_AL_CUT10_comp
         COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_numerical.py
         ${CMAKE_SOURCE_DIR}/tests/1ao6/ENE_1_WANG_AL_CUT10.ref
         Testing/output 
         ${RTOL} ${ATOL})

# Tests on the full potential
## TODO ## add_test(NAME NMA_ene_full_AMOEBA
## TODO ##          COMMAND  bin/C_test_potential_xyz
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
## TODO ##          amoeba09.prm
## TODO ##          Testing/output)
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_eb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ev
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EV
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_em
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EM
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ep
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EP
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_eub
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EUB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ea
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_eba
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_eopb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EOPB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ept
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EPT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_et
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ET
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ett
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_eat
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EAT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_ebt
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME NMA_ene_full_AMOEBA_comp_etot
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETOT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## 
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA
## TODO ##          COMMAND  bin/C_test_potential_xyz
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz
## TODO ##          amoebabio18.prm
## TODO ##          Testing/output)
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_eb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ev
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EV
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_em
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EM
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ep
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EP
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_eub
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EUB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ea
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_eba
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_eopb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EOPB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ept
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EPT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_et
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ET
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ett
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_eat
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EAT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_ebt
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1crn_ene_full_AMOEBA_comp_etot
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETOT
## TODO ##          2.0e-5 3.0e-4)
## TODO ## 
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA
## TODO ##          COMMAND  bin/C_test_potential_xyz
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz
## TODO ##          amoebabio18.prm
## TODO ##          Testing/output)
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_eb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ev
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EV
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_em
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EM
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ep
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EP
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_eub
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EUB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ea
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_eba
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_eopb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EOPB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ept
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EPT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_et
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ET
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ett
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_eat
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EAT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_ebt
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 1ubq_ene_full_AMOEBA_comp_etot
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETOT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## 
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA
## TODO ##          COMMAND  bin/C_test_potential_xyz
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/input.xyz
## TODO ##          amoebabio18.prm
## TODO ##          Testing/output)
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_eb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ev
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EV
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_em
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EM
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ep
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EP
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_eub
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EUB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ea
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EA
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_eba
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBA
## TODO ##          1e-4 1e-5)
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_eopb
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EOPB
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ept
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EPT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_et
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ET
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ett
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_eat
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EAT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_ebt
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          EBT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## add_test(NAME 2d94_ene_full_AMOEBA_comp_etot
## TODO ##          COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##          ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##          Testing/output 
## TODO ##          ETOT
## TODO ##          ${RTOL} ${ATOL})
## TODO ## 
## TODO ## # Tests on the full potential (after reloading from HDF5)
## TODO ## if (HDF5_WORKS)
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_hdf5_dump
## TODO ##              COMMAND  bin/ommp_pp
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/input.xyz
## TODO ##              amoeba09.prm
## TODO ##              Testing/NMA.hdf5)
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_hdf5_load
## TODO ##              COMMAND  bin/C_test_potential_hdf5
## TODO ##              Testing/NMA.hdf5
## TODO ##              Testing/output_hdf5)
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_eb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5
## TODO ##              EB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ev_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EV
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_em_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EM
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ep_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EP
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_eub_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EUB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ea_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_eba_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_eopb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EOPB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ept_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EPT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_et_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ET
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ett_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_eat_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EAT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_ebt_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME NMA_ene_full_AMOEBA_comp_etot_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/N-methylacetamide/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETOT
## TODO ##              ${RTOL} ${ATOL})
## TODO ## 
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_hdf5_dump
## TODO ##              COMMAND  bin/ommp_pp
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/input.xyz
## TODO ##              amoebabio18.prm
## TODO ##              Testing/1crn.hdf5)
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_hdf5_load
## TODO ##              COMMAND  bin/C_test_potential_hdf5
## TODO ##              Testing/1crn.hdf5
## TODO ##              Testing/output_hdf5)
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_eb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ev_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EV
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_em_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EM
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ep_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EP
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_eub_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EUB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ea_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_eba_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_eopb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EOPB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ept_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EPT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_et_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ET
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ett_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_eat_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EAT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_ebt_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1crn_ene_full_AMOEBA_comp_etot_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1crn/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETOT
## TODO ##              2.0e-5 3.0e-4)
## TODO ## 
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_hdf5_dump
## TODO ##              COMMAND  bin/ommp_pp
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/input.xyz
## TODO ##              amoebabio18.prm
## TODO ##              Testing/1ubq.hdf5)
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_hdf5_load
## TODO ##              COMMAND  bin/C_test_potential_hdf5
## TODO ##              Testing/1ubq.hdf5
## TODO ##              Testing/output_hdf5)
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ev_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EV
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_em_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EM
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ep_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EP
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eub_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EUB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ea_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eba_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eopb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EOPB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ept_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EPT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_et_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ET
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ett_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_eat_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EAT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_ebt_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 1ubq_ene_full_AMOEBA_comp_etot_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/1ubq/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETOT
## TODO ##              ${RTOL} ${ATOL})
## TODO ## 
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_hdf5_dump
## TODO ##              COMMAND  bin/ommp_pp
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/input.xyz
## TODO ##              amoebabio18.prm
## TODO ##              Testing/2d94.hdf5)
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_hdf5_load
## TODO ##              COMMAND  bin/C_test_potential_hdf5
## TODO ##              Testing/2d94.hdf5
## TODO ##              Testing/output_hdf5)
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_eb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ev_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EV
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_em_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EM
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ep_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EP
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_eub_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EUB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ea_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EA
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_eba_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBA
## TODO ##              1e-4 1e-5)
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_eopb_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EOPB
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ept_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EPT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_et_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ET
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ett_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_eat_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EAT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_ebt_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              EBT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##     add_test(NAME 2d94_ene_full_AMOEBA_comp_etot_hdf5
## TODO ##              COMMAND python3 ${CMAKE_SOURCE_DIR}/tests/compare_full_pot.py
## TODO ##              ${CMAKE_SOURCE_DIR}/tests/2d94/FULL_POTENTIAL.ref
## TODO ##              Testing/output_hdf5 
## TODO ##              ETOT
## TODO ##              ${RTOL} ${ATOL})
## TODO ##endif ()
