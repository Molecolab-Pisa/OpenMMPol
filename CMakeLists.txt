cmake_minimum_required(VERSION 3.17)
project(Open-MMPol)
enable_language(Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_options($<$<CONFIG:RELEASE>:-O3>
                    $<$<CONFIG:DEBUG>:-Og>
                    $<$<CONFIG:DEBUG>:-g>)

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
            src/elstat.f90
            src/energy.f90
            src/polar.f90
            src/polarization.f90
            src/rotate_multipoles.f90
            src/solvers.f90
            src/utilities.f90)

# LAPACK
find_package(LAPACK REQUIRED COMPONENTS Fortran)
target_link_libraries(openmmpol ${LAPACK_LIBRARIES})
# cmake 3.22 exposes BLA_SIZEOF_INTEGER

# HDF5
find_package(HDF5 COMPONENTS Fortran)
if (HDF5_FOUND)
    add_compile_definitions(USE_HDF5)
    include_directories( SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(openmmpol ${HDF5_Fortran_LIBRARIES})
endif(HDF5_FOUND)

target_compile_options(openmmpol PRIVATE -cpp)
# cmake 3.18 is needed for this
#set_source_files_properties(src/mod_interface.f90 src/mod_io.f90 src/mod_memory.f90 PROPERTIES Fortran_PREPROCESS ON)

target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2003)
target_compile_options(openmmpol PRIVATE -fPIC)
target_compile_options(openmmpol PRIVATE -fall-intrinsics)
target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-fbacktrace>)

include_directories(include)
add_executable(test_c src/test.c)
target_link_libraries(test_c openmmpol)
target_compile_options(test_c PRIVATE -Wall -Wextra -pedantic)
set_property(TARGET test_c PROPERTY C_STANDARD 99)

