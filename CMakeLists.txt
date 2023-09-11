cmake_minimum_required(VERSION 3.17)
project(openMMPol)
enable_language(Fortran)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_options($<$<CONFIG:RELEASE>:-O3>
                    $<$<CONFIG:DEBUG>:-g>
                    $<$<CONFIG:COVERAGE>:-O0>)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    add_compile_options($<$<CONFIG:RELEASE>:>
                        $<$<CONFIG:DEBUG>:-Og>
                        $<$<CONFIG:COVERAGE>:-pg>
                        $<$<CONFIG:COVERAGE>:-fprofile-arcs>
                        $<$<CONFIG:COVERAGE>:-ftest-coverage>)
    add_link_options($<$<CONFIG:COVERAGE>:-coverage>)
endif()

include_directories(include) # Directory for .h files

# Create the define version.h needed for tracking version
# of the library.
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} -C ${CMAKE_SOURCE_DIR} rev-parse OUTPUT_QUIET ERROR_QUIET RESULT_VARIABLE CMAKE_SOURCE_DIR_IS_GIT)

if (CMAKE_SOURCE_DIR_IS_GIT EQUAL 0)
    message(NOTICE "Version read from git repository information.")
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --dirty
                    COMMAND sed "s/-/+r/;s/-/./g"
                    OUTPUT_VARIABLE OMMP_VERSION)
    string(STRIP ${OMMP_VERSION} OMMP_VERSION)
else()
    if(NOT EXISTS ${CMAKE_SOURCE_DIR}/include/version.h)
        set(OMMP_VERSION "unknown-nogit")
        message(DEPRECATION "Your source folder is not a git repository, so an accurate tracking of version cannot be performed. Version name is \"${OMMP_VERSION}\". We strongly recommend to install OpenMMPol from a git repository to have a proper version identifier needed for reproducing results.")
    else()
        message(NOTICE "Version read from ${CMAKE_SOURCE_DIR}/include/version.h")
        execute_process(COMMAND echo \#include \"version.h\"\n_OMMP_VERSION
                        COMMAND ${CMAKE_C_COMPILER} -E -P -I${CMAKE_SOURCE_DIR}/include -
                        COMMAND tr -d \"
                        OUTPUT_VARIABLE OMMP_VERSION)
    endif()
endif()

file(WRITE ${CMAKE_BINARY_DIR}/include/version.h "#define _OMMP_VERSION \"${OMMP_VERSION}\"")
configure_file(${CMAKE_SOURCE_DIR}/include/_openmmpol_const.h ${CMAKE_BINARY_DIR}/include/openmmpol_const.h)
file(COPY ${CMAKE_SOURCE_DIR}/include/openmmpol.h DESTINATION ${CMAKE_BINARY_DIR}/include/)
include_directories(${CMAKE_BINARY_DIR}/include) # Directory for .h files
message(NOTICE "Version string ${OMMP_VERSION}")

set(LIBRARY_SRC
    src/mod_adjacency_mat.f90
    src/mod_bonded.f90
    src/mod_constants.f90
    src/mod_c_interface.f90
    src/mod_electrostatics.f90
    src/mod_geomgrad.f90
    src/mod_inputloader.f90
    src/mod_interface.f90
    src/mod_io.f90
    src/mod_iohdf5.f90
    src/mod_jacobian_mat.f90
    src/mod_link_atom.f90
    src/mod_memory.f90
    src/mod_mmpol.f90
    src/mod_nonbonded.f90
    src/mod_polarization.f90
    src/mod_prm.f90
    src/mod_qm_helper.f90
    src/mod_solvers.f90
    src/mod_topology.f90
    src/mod_utils.f90
    src/rotate_multipoles.f90
    src/smartinput.c)


add_library(openmmpol SHARED ${LIBRARY_SRC}) 

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

# cJSON
find_package(cJSON REQUIRED)
if(cJSON_FOUND)
    include_directories( SYSTEM ${CJSON_INCLUDE_DIR})
    target_link_libraries(openmmpol ${CJSON_LIBRARIES})
endif()

# openssl
find_package(OpenSSL REQUIRED)
if(OpenSSL_FOUND)
    include_directories( SYSTEM ${OPENSSL_INCLUDE_DIR})
    target_link_libraries(openmmpol ${OPENSSL_LIBRARIES})
endif()

# Enable preprocessor for Fortran
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
# cmake 3.18 is needed for this
#set_source_files_properties(src/mod_interface.f90 src/mod_io.f90 src/mod_memory.f90 PROPERTIES Fortran_PREPROCESS ON)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008ts -fno-realloc-lhs -fall-intrinsics")
    target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic) # Very strict check on code
    #target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2003 -fno-realloc-lhs) # Very strict check on code
    target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-g -Wconversion>)
    if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Warray-temporaries -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ffree-line-length-0 -fcheck=all")
    endif()
    target_link_libraries(openmmpol stdc++)
endif()
target_compile_options(openmmpol PRIVATE -fPIC) # Generate position independent code, for library

set(TEST_PROGRAMS
    test_timing
    test_SI_init
    test_SI_potential
    test_SI_geomgrad
    test_SI_geomgrad_num)

set(UTILITY_PROGRAMS 
    ommp_xyz2mmp)

if (HDF5_WORKS)
    set(UTILITY_PROGRAMS
        ${UTILITY_PROGRAMS}
        ommp_pp)

    set(TEST_PROGRAMS
        ${TEST_PROGRAMS}
        test_init_hdf5
        test_potential_hdf5)
endif ()

if(NOT PYTHON)
    set(UTILITY_EXE)
    foreach( exe ${UTILITY_PROGRAMS} )
        add_executable(${exe} src/${exe}.c)
        target_compile_options(${exe} PRIVATE -Wall -Wextra -pedantic)
        set_property(TARGET ${exe} PROPERTY C_STANDARD 99)
        target_link_libraries(${exe} openmmpol)
        set(UTILITY_EXE
            ${UTILITY_EXE}
            ${exe})
    endforeach()
    add_custom_target(utility_programs DEPENDS ${UTILITY_PROGRAMS})

    set(TEST_LANG
        C
        F03)
    foreach( lang ${TEST_LANG} )
        set(TEST_EXE)
        foreach( exe ${TEST_PROGRAMS} )
            if (${lang} STREQUAL "C")
                if(NOT EXISTS ${CMAKE_SOURCE_DIR}/src/test_programs/${lang}/${exe}.c)
                    continue ()
                endif ()

                add_executable(${lang}_${exe} src/test_programs/${lang}/${exe}.c)
                target_compile_options(${lang}_${exe} PRIVATE -Wall -Wextra -pedantic)
                set_property(TARGET ${lang}_${exe} PROPERTY C_STANDARD 99)
            endif()
            if (${lang} STREQUAL "F03") 
                if(NOT EXISTS ${CMAKE_SOURCE_DIR}/src/test_programs/${lang}/${exe}.f90) 
                    continue ()
                endif ()

                add_executable(${lang}_${exe} src/test_programs/${lang}/${exe}.f90)
                if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
                    target_compile_options(${lang}_${exe} PRIVATE -Wall -Wextra -pedantic -std=f2003)
                endif ()
            endif()

            target_link_libraries(${lang}_${exe} openmmpol)
            set(TEST_EXE
                ${TEST_EXE}
                ${lang}_${exe})
        endforeach()
        add_custom_target(${lang}_test_programs DEPENDS ${TEST_EXE})
    endforeach()
endif()

# Python
option(PYTHON "Build python module (Generally you want to use ./setup.py build_ext for this)" OFF)
if (PYTHON)
    #python
    find_package(pybind11 CONFIG REQUIRED)
    add_library(openmmpol_static STATIC
                ${LIBRARY_SRC})
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
        target_compile_options(openmmpol_static PRIVATE -Wall -Wextra -pedantic -std=f2008ts -fno-realloc-lhs)
        target_compile_options(openmmpol_static PRIVATE -fall-intrinsics) # This is only needed for a sizeof function
    endif()
    target_compile_options(openmmpol_static PRIVATE -fPIC -cpp)
    target_link_libraries(openmmpol_static OpenMP::OpenMP_Fortran)
    target_link_libraries(openmmpol_static ${LAPACK_LIBRARIES})
    if (HDF5_WORKS)
        target_link_libraries(openmmpol_static ${HDF5_Fortran_LIBRARIES})
        add_definitions(-DUSE_HDF5)
    endif()
    target_link_libraries(openmmpol_static ${CJSON_LIBRARIES})
    target_link_libraries(openmmpol_static ${OPENSSL_LIBRARIES})

    pybind11_add_module(pyopenmmpol MODULE src/pyommp_interface.cpp)
    target_compile_definitions(pyopenmmpol PRIVATE VERSION_INFO=${OMMP_VERSION})
    target_link_libraries(pyopenmmpol PRIVATE openmmpol_static)

    set_target_properties(pyopenmmpol PROPERTIES LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    set_target_properties(pyopenmmpol PROPERTIES OUTPUT_NAME pyopenmmpol)
endif()

include(GNUInstallDirs)
install(TARGETS C_test_timing openmmpol ${UTILITY_EXE}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/openMMPol-${OMMP_VERSION}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/openMMPol-${OMMP_VERSION})
install(TARGETS C_test_timing ${UTILITY_EXE})
install(FILES ${CMAKE_BINARY_DIR}/include/openmmpol.h ${CMAKE_BINARY_DIR}/mod/ommp_interface.mod
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/openMMPol-${OMMP_VERSION})

configure_file(${CMAKE_SOURCE_DIR}/cmake/openmmpol-config.cmake.in 
               ${CMAKE_BINARY_DIR}/cmake/openmmpol-config.cmake @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/openmmpol-config-version.cmake.in 
               ${CMAKE_BINARY_DIR}/cmake/openmmpol-config-version.cmake @ONLY)

install(FILES ${CMAKE_BINARY_DIR}/cmake/openmmpol-config.cmake
              ${CMAKE_BINARY_DIR}/cmake/openmmpol-config-version.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/openMMPol-${OMMP_VERSION})
# Tests
include(tests/CMakeLists.txt)
