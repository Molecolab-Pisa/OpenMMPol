#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2008ts -fno-realloc-lhs -fall-intrinsics -Wall -Wextra -pedantic"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "-g ${Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "-g -Wall -std=f2008ts -fbounds-check "
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

# if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_LESS 8.0)
#     set(OpenMP_Fortran_FLAGS "-fopenmp")
#     set(OpenMP_Fortran_LIB_NAMES "${OpenMP_C_LIB_NAMES}")
# endif()
