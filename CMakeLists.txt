cmake_minimum_required(VERSION 3.17)
project(Open-MMPol)
enable_language(Fortran)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_options($<$<CONFIG:RELEASE>:-O3>
                    $<$<CONFIG:DEBUG>:-g>)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    add_compile_options($<$<CONFIG:DEBUG>:-Og>)
endif()

#set(CMAKE_VERBOSE_MAKEFILE ON)

add_library(openmmpol SHARED 
            src/mod_adjacency_mat.f90
            #src/mod_bonded.f90
            src/mod_constants.f90
            src/mod_electrostatics.f90
            #src/mod_inputloader.f90
            #src/mod_interface.f90
            src/mod_io.f90
            #src/mod_iohdf5.f90
            src/mod_memory.f90
            src/mod_mmpol.f90
            #src/mod_nonbonded.f90
            #src/mod_polarization.f90
            #src/mod_prm.f90
            #src/mod_solvers.f90
            src/mod_topology.f90
            src/mod_utils.f90
            src/rotate_multipoles.f90)

# LAPACK
find_package(LAPACK REQUIRED COMPONENTS Fortran)
target_link_libraries(openmmpol ${LAPACK_LIBRARIES})
# cmake 3.22 exposes BLA_SIZEOF_INTEGER

#OpenMP
find_package(OpenMP REQUIRED COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
    include_directories(SYSTEM ${OpenMP_Fortran_INCLUDE_DIRS})
    target_link_libraries(openmmpol OpenMP::OpenMP_Fortran)
endif()

# HDF5
find_package(HDF5 COMPONENTS Fortran)
if (HDF5_FOUND)
    include_directories( SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(openmmpol ${HDF5_Fortran_LIBRARIES})
    set(CMAKE_REQUIRED_INCLUDES ${HDF5_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${HDF5_Fortran_LIBRARIES})
    include(CheckFortranSourceCompiles)
    check_fortran_source_compiles("program test_hdf5; use hdf5; integer(hid_t) :: test; integer(kind=4) :: eflag; call h5open_f(eflag); end"
                                  HDF5_WORKS SRC_EXT .f90)
    if (HDF5_WORKS)
        add_definitions(-DUSE_HDF5)
    endif() 
endif()

# Enable preprocessor for Fortran
target_compile_options(openmmpol PRIVATE -cpp)
# cmake 3.18 is needed for this
#set_source_files_properties(src/mod_interface.f90 src/mod_io.f90 src/mod_memory.f90 PROPERTIES Fortran_PREPROCESS ON)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2008ts -fno-realloc-lhs) # Very strict check on code
    target_compile_options(openmmpol PRIVATE -fall-intrinsics) # This is only needed for a sizeof function
    target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-g -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan>)
endif()
target_compile_options(openmmpol PRIVATE -fPIC) # Generate position independent code, for library

#include_directories(include) # Directory for .h files
#
#add_executable(C_test_init_xyz src/test_programs/test_init_xyz.c)
#target_link_libraries(C_test_init_xyz openmmpol)
#target_compile_options(C_test_init_xyz PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_init_xyz PROPERTY C_STANDARD 99)
#
#add_executable(C_test_init src/test_programs/test_init.c)
#target_link_libraries(C_test_init openmmpol)
#target_compile_options(C_test_init PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_init PROPERTY C_STANDARD 99)
#
#add_executable(C_test_reinit src/test_programs/test_reinit.c)
#target_link_libraries(C_test_reinit openmmpol)
#target_compile_options(C_test_reinit PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_reinit PROPERTY C_STANDARD 99)
#
#add_executable(C_test_ipd src/test_programs/test_ipd.c)
#target_link_libraries(C_test_ipd openmmpol)
#target_compile_options(C_test_ipd PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_ipd PROPERTY C_STANDARD 99)
#
#add_executable(C_test_potential_xyz src/test_programs/test_potential_xyz.c)
#target_link_libraries(C_test_potential_xyz openmmpol)
#target_compile_options(C_test_potential_xyz PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_potential_xyz PROPERTY C_STANDARD 99)
#
#add_executable(C_test_energy src/test_programs/test_energy.c)
#target_link_libraries(C_test_energy openmmpol)
#target_compile_options(C_test_energy PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET C_test_energy PROPERTY C_STANDARD 99)
#
#add_executable(ommp_xyz2mmp src/ommp_xyz2mmp.c)
#target_link_libraries(ommp_xyz2mmp openmmpol)
#target_compile_options(ommp_xyz2mmp PRIVATE -Wall -Wextra -pedantic)
#set_property(TARGET ommp_xyz2mmp PROPERTY C_STANDARD 99)
#
#if (HDF5_WORKS)
#    add_executable(ommp_pp src/ommp_compile.c)
#    target_link_libraries(ommp_pp openmmpol)
#    target_compile_options(ommp_pp PRIVATE -Wall -Wextra -pedantic)
#    set_property(TARGET ommp_pp PROPERTY C_STANDARD 99)
#    
#    add_executable(C_test_init_hdf5 src/test_programs/test_init_hdf5.c)
#    target_link_libraries(C_test_init_hdf5 openmmpol)
#    target_compile_options(C_test_init_hdf5 PRIVATE -Wall -Wextra -pedantic)
#    set_property(TARGET C_test_init_hdf5 PROPERTY C_STANDARD 99)
#    
#    add_executable(C_test_potential_hdf5 src/test_programs/test_potential_hdf5.c)
#    target_link_libraries(C_test_potential_hdf5 openmmpol)
#    target_compile_options(C_test_potential_hdf5 PRIVATE -Wall -Wextra -pedantic)
#    set_property(TARGET C_test_potential_hdf5 PROPERTY C_STANDARD 99)
#endif ()
#
#
#add_executable(F_test_init src/test_programs/test_init.F90)
#target_link_libraries(F_test_init openmmpol)
#if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#    target_compile_options(F_test_init PRIVATE -Wall -Wextra -pedantic -std=f2003)
#endif()
#
#add_executable(F_test_ipd src/test_programs/test_ipd.F90)
#target_link_libraries(F_test_ipd openmmpol)
#if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#    target_compile_options(F_test_ipd PRIVATE -Wall -Wextra -pedantic -std=f2003)
#endif()
#
#include(tests/CMakeLists.txt)

add_executable(F_test_tmp src/test_programs/test_tmp.F90)
target_link_libraries(F_test_tmp openmmpol)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    target_compile_options(F_test_tmp PRIVATE -Wall -Wextra -pedantic -std=f2003)
endif()
