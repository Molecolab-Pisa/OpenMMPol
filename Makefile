#
# makefile for open_mmpol
#
RunF77 = gfortran
FFLAGS = -cpp -O3 -fopenmp -march=native -g -fbacktrace
BLAS   = -lblas -llapack
#
MODS   = precision.o mmpol.o solvers.o multipoles_functions.o polarization_functions.o
OBJS   = coulomb_kernel.o electrostatics.o main_amoeba.o mmpol_init.o mmpol_process.o rotate_multipoles.o utilities.o
#
all:    $(MODS) $(OBJS)
	$(RunF77) $(FFLAGS) -o main.exe $(MODS) $(OBJS) $(BLAS)

#
%.o: %.f
	$(RunF77) $(FFLAGS) -c $*.f
%.o: %.f90
	$(RunF77) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
