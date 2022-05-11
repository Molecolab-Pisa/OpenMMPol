# makefile for open_mmpol

DEBUG = YES
#--------------------------------------#
OBJ_DIR = ./obj
BIN_DIR = ./bin
LIB_DIR = ./lib
MOD_DIR = ./mod
SRC_DIR = ./src
PYT_DIR = ./python/pyopenmmp/pyopenmmp

FC = gfortran

ifeq ($(DEBUG), YES)
	FFLAGS = -Wall -Wextra -pedantic -std=f2003 -fcheck=all -g -Og -fbacktrace -fPIC
else
	FFLAGS = -Wall -Wextra -pedantic -std=f2003 -O3 -fPIC
endif

CPPFLAGS = -cpp
LDLIBS = -lblas -llapack -lgfortran

OBJS   = coulomb_kernel.o \
	 electrostatics.o \
	 mmpol.o \
	 mmpol_init.o \
	 mmpol_process.o \
	 elstat.o \
	 energy.o \
         polar.o \
	 polarization.o \
	 precision.o \
	 rotate_multipoles.o \
	 solvers.o \
	 utilities.o

OBJS := $(addprefix $(OBJ_DIR)/, $(OBJS))

PY_SUFFIX = $(shell python3-config --extension-suffix)

all: libraries binaries python

python: $(PYT_DIR)/pymmpol.so

libraries: $(LIB_DIR)/mmpolmodules.a $(LIB_DIR)/libopenmmpol.so

binaries: $(BIN_DIR)/main.exe $(BIN_DIR)/main_amoeba.exe

$(BIN_DIR)/%.exe: $(SRC_DIR)/%.F90 $(BIN_DIR) $(LIB_DIR)/libopenmmpol.so
	$(FC) $(CPPFLAGS) $(FFLAGS) $(LDLIBS) -L$(LIB_DIR) -I$(MOD_DIR) $< -lopenmmpol -o $@

$(LIB_DIR)/mmpolmodules.a: $(OBJS) $(LIB_DIR)
	ar crs $@ $(OBJS)

$(LIB_DIR)/libopenmmpol.so: $(OBJS) $(LIB_DIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) $(LDLIBS) -shared -I$(MOD_DIR) -o $@ $(OBJS)

$(PYT_DIR)/pymmpol.so: $(PYT_DIR)/pymmpol$(PY_SUFFIX)
	ln -s $(shell realpath $<) $@

$(PYT_DIR)/pymmpol$(PY_SUFFIX): $(OBJS) $(LIB_DIR)/mmpolmodules.a
	echo "{'real': {'rp': 'double'}, 'integer': {'ip': 'long'}}" > .f2py_f2cmap
	f2py3 -c -lblas -llapack -m pymmpol src/precision.f90 \
					    src/wrapper.f90 \
		                            src/mmpol.f90 \
					    src/mmpol_init.f90 \
					    lib/mmpolmodules.a \
					    skip: r_alloc1 r_alloc2 r_alloc3 \
					          i_alloc1 i_alloc2 i_alloc3 \
						  r_free1 r_free2 r_free3 \
						  i_free1 i_free2 i_free3 : \
					    --fcompiler=gnu95 \
					    --f90flags="$(FFLAGS)"
	mv pymmpol$(PY_SUFFIX) $(PYT_DIR)
	rm .f2py_f2cmap

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
$(MOD_DIR):
	@mkdir -p $(MOD_DIR)
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)
$(LIB_DIR):
	@mkdir -p $(LIB_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(MOD_DIR) $(LIB_DIR)
	rm -rf $(PYT_DIR)/*.so

# Explicit dependencies
$(OBJ_DIR)/coulomb_kernel.o: $(OBJ_DIR)/mmpol.o
$(OBJ_DIR)/elstat.o:
$(OBJ_DIR)/electrostatics.o: $(OBJ_DIR)/elstat.o $(OBJ_DIR)/mmpol.o $(OBJ_DIR)/precision.o
$(OBJ_DIR)/energy.o: $(OBJ_DIR)/mmpol.o
$(OBJ_DIR)/mmpol.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/mmpol_init.o: $(OBJ_DIR)/mmpol.o
$(OBJ_DIR)/mmpol_process.o: $(OBJ_DIR)/mmpol.o
$(OBJ_DIR)/multipoles_functions.o: $(OBJ_DIR)/elstat.o 
$(OBJ_DIR)/polar.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/polarization.o: $(OBJ_DIR)/solvers.o
$(OBJ_DIR)/precision.o:
$(OBJ_DIR)/rotate_multipoles.o: $(OBJ_DIR)/mmpol.o 
$(OBJ_DIR)/solvers.o: $(OBJ_DIR)/precision.o
$(OBJ_DIR)/utilities.o: $(OBJ_DIR)/mmpol.o
