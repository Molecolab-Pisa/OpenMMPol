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

# set(CMAKE_VERBOSE_MAKEFILE ON)

add_library(openmmpol SHARED 
            src/mod_adjacency_mat.f90
            src/mod_bicubic_interp.f90
            src/mod_bonded.f90
            src/mod_constants.f90
            src/mod_electrostatics.f90
            src/mod_inputloader.f90
            src/mod_interface.f90
            src/mod_io.f90
            src/mod_memory.f90
            src/mod_mmpol.f90
            src/mod_nonbonded.f90
            src/mod_polarization.f90
            src/mod_prm.f90
            src/mod_solvers.f90
            src/mod_utils.f90
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
target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-g -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan>)

include_directories(include) # Directory for .h files

add_executable(C_test_init_xyz src/test_init_xyz.c)
target_link_libraries(C_test_init_xyz openmmpol)
target_compile_options(C_test_init_xyz PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_init_xyz PROPERTY C_STANDARD 99)

add_executable(C_test_init src/test_init.c)
target_link_libraries(C_test_init openmmpol)
target_compile_options(C_test_init PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_init PROPERTY C_STANDARD 99)

add_executable(C_test_ipd src/test_ipd.c)
target_link_libraries(C_test_ipd openmmpol)
target_compile_options(C_test_ipd PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_ipd PROPERTY C_STANDARD 99)

add_executable(C_test_potential_xyz src/test_potential_xyz.c)
target_link_libraries(C_test_potential_xyz openmmpol)
target_compile_options(C_test_potential_xyz PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_potential_xyz PROPERTY C_STANDARD 99)

add_executable(C_test_energy src/test_energy.c)
target_link_libraries(C_test_energy openmmpol)
target_compile_options(C_test_energy PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET C_test_energy PROPERTY C_STANDARD 99)

add_executable(F_test_init src/test_init.F90)
target_link_libraries(F_test_init openmmpol)
target_compile_options(F_test_init PRIVATE -Wall -Wextra -pedantic -std=f2003)

include(tests/CMakeLists.txt)
