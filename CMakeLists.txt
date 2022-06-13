cmake_minimum_required(VERSION 3.17)
project(Open-MMPol)
enable_language(Fortran)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_options($<$<CONFIG:RELEASE>:-O3>
                    $<$<CONFIG:DEBUG>:-Og>
                    $<$<CONFIG:DEBUG>:-g>)

#set(CMAKE_VERBOSE_MAKEFILE ON)

add_library(openmmpol SHARED 
            src/coulomb_kernel.f90 
            src/electrostatics.f90
            src/mod_adjacency_mat.f90
            src/mod_constants.f90
            src/mod_inputloader.f90
            src/mod_interface.f90
            src/mod_io.f90
            src/mod_memory.f90
            src/mod_mmpol.f90
            src/mod_polarization.f90
            src/mod_solvers.f90
            src/elstat.f90
            src/energy.f90
            src/rotate_multipoles.f90)

# LAPACK
find_package(LAPACK REQUIRED COMPONENTS Fortran)
target_link_libraries(openmmpol ${LAPACK_LIBRARIES})
# cmake 3.22 exposes BLA_SIZEOF_INTEGER

# HDF5
find_package(HDF5 COMPONENTS Fortran)
if (HDF5_FOUND)
    add_definitions(-DUSE_HDF5)
    include_directories( SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(openmmpol ${HDF5_Fortran_LIBRARIES})
endif(HDF5_FOUND)

# Enable preprocessor for Fortran
target_compile_options(openmmpol PRIVATE -cpp)
# cmake 3.18 is needed for this
#set_source_files_properties(src/mod_interface.f90 src/mod_io.f90 src/mod_memory.f90 PROPERTIES Fortran_PREPROCESS ON)

target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2003) # Very strict check on code
target_compile_options(openmmpol PRIVATE -fPIC) # Generate position independent code, for library
target_compile_options(openmmpol PRIVATE -fall-intrinsics) # This is only needed for a sizeof function
target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-fbacktrace>)

include_directories(include) # Directory for .h files

add_executable(C_test_init src/test_init.c)
target_link_libraries(C_test_init openmmpol)
target_compile_options(C_test_init PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_init PROPERTY C_STANDARD 99)

add_executable(C_test_ipd src/test_ipd.c)
target_link_libraries(C_test_ipd openmmpol)
target_compile_options(C_test_ipd PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_ipd PROPERTY C_STANDARD 99)

add_executable(C_test_energy src/test_energy.c)
target_link_libraries(C_test_energy openmmpol)
target_compile_options(C_test_energy PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_energy PROPERTY C_STANDARD 99)

add_executable(F_test_init src/test_init.F90)
target_link_libraries(F_test_init openmmpol)
target_compile_options(F_test_init PRIVATE -Wall -Wextra -pedantic -std=f2003)

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
