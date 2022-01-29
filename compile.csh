#!/bin/tcsh

# Temporary compile script for f2py

set MODS = (precision.o mmpol.o solvers.o multipoles_functions.o polarization_functions.o)
set OBJS = (coulomb_kernel.o electrostatics.o main_amoeba.o mmpol_init.o mmpol_process.o rotate_multipoles.o utilities.o polarization.o energy.o)

set skip = (r_alloc1 r_alloc2 r_alloc3 i_alloc1 i_alloc2 i_alloc3 r_free1 r_free2 r_free3 i_free1 i_free2 i_free3)

# Compile -fPIC
foreach file ($MODS $OBJS)
  gfortran -cpp -O3 -fopenmp -march=native -g -fbacktrace -fPIC -fdefault-integer-8 -c ${file:r}.f90 
end


# Make .a library
ar crs mmpolmodules.a $MODS $OBJS

# Force precision
echo "{'real': {'rp': 'double'}, 'integer': {'ip': 'long'}}" > .f2py_f2cmap

f2py3 -c -lblas -llapack -m pymmpol wrapper.f90 mmpol.f90 mmpol_init.f90 mmpolmodules.a skip: $skip  : --fcompiler=gnu95 --f90flags="-cpp -O3 -fdefault-integer-8 -fopenmp -march=native -g -fbacktrace -fPIC" 



