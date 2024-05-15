#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g -Wall -std=f2018 -fbounds-check"
  CACHE STRING "Fortran compiler flags for Debug build")

set(Fortran_FLAGS_COVERAGE "-O0 -g --coverage")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "-O2 -funroll-all-loops"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_RELWITDEBINFO "-g ${C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for RelWithDebInfo build")

set(C_FLAGS_DEBUG "-g -Wall -pedantic -fbounds-check"
  CACHE STRING "C compiler flags for Debug build")

set(C_FLAGS_COVERAGE "-O0 -g --coverage")

# TODO: ugly, refactor
# if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
#   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2018 -fno-realloc-lhs -fall-intrinsics")
#   target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic) # Very strict check on code
#   #target_compile_options(openmmpol PRIVATE -Wall -Wextra -pedantic -std=f2003 -fno-realloc-lhs) # Very strict check on code
#   target_compile_options(openmmpol PRIVATE $<$<CONFIG:DEBUG>:-g -Wconversion>)
#   if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
#       set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Warray-temporaries -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan -ffree-line-length-0 -fcheck=all")
#   endif()
#   # target_link_libraries(openmmpol stdc++)
#   set (lib-deps 
#       "${lib-deps}"
#       "stdc++")
# endif()