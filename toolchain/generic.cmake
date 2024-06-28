#
# Fortran compiler settings
#
set(Fortran_FLAGS "${CMAKE_Fortran_FLAGS}"
  CACHE STRING "Build type independent Fortran compiler flags")

set(Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_RELWITHDEBINFO "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
  CACHE STRING "Fortran compiler flags for Release build")

set(Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG}"
  CACHE STRING "Fortran compiler flags for Debug build")

# Use intrinsic Fortran 2008 erf/erfc functions
set(INTERNAL_ERFC CACHE BOOL 0)

#
# C compiler settings
#
set(C_FLAGS "${CMAKE_C_FLAGS}"
  CACHE STRING "Build type independent C compiler flags")

set(C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}"
  CACHE STRING  "C compiler flags for Release build")

set(C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}"
  CACHE STRING "C compiler flags for Debug build")