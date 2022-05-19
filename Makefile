# makefile for open_mmpol

DEBUG = YES
USE_INT64 = NO
USE_HDF5 = YES
#--------------------------------------#
DOC_DIR = ./doc
OBJ_DIR = ./obj
BIN_DIR = ./bin
LIB_DIR = ./lib
MOD_DIR = ./mod
SRC_DIR = ./src
INCLUDE_DIR = ./include

FC = gfortran
CC = gcc

ifeq ($(DEBUG), YES)
	FFLAGS = -Wall -Wextra -pedantic -std=f2003 -fcheck=all -fall-intrinsics -g -Og -fbacktrace -fPIC
else
	FFLAGS = -Wall -Wextra -pedantic -std=f2003 -O3 -fPIC
endif

CPPFLAGS = -cpp
CFLAGS = -Wall -pedantic -g -Og -std=c99
LDLIBS = -lblas -llapack -lgfortran
DOC_CPPFLAGS = 

ifeq ($(USE_INT64), YES)
    CPPFLAGS += -DUSE_I8 
	DOC_CPPFLAGS += -m USE_I8
endif

ifeq ($(USE_HDF5), YES)
    CPPFLAGS += -DUSE_HDF5 
	DOC_CPPFLAGS += -m USE_HDF5
	LDLIBS += -lhdf5_fortran
	FFLAGS += -I/usr/include
endif


OBJS   = coulomb_kernel.o \
	     electrostatics.o \
	     mmpol_init.o \
         mod_constants.o \
         mod_interface.o \
         mod_io.o \
         mod_memory.o \
	     mod_mmpol.o \
	     elstat.o \
	     energy.o \
         polar.o \
	     polarization.o \
	     rotate_multipoles.o \
	     solvers.o \
	     utilities.o

OBJS := $(addprefix $(OBJ_DIR)/, $(OBJS))

all: libraries binaries # python

python: $(PYT_DIR)/pymmpol.so

libraries: $(LIB_DIR)/mmpolmodules.a $(LIB_DIR)/libopenmmpol.so

binaries: $(BIN_DIR)/test_init.exe $(BIN_DIR)/test_amoeba.exe $(BIN_DIR)/test_c.out

doc: openMMPol.md
	ford openMMPol.md $(DOC_CPPFLAGS)

$(BIN_DIR)/%.exe: $(SRC_DIR)/%.F03 $(BIN_DIR) $(LIB_DIR)/libopenmmpol.so
	$(FC) $(CPPFLAGS) $(FFLAGS) $(LDLIBS) -L$(LIB_DIR) -I$(MOD_DIR) $< -lopenmmpol -o $@

$(BIN_DIR)/%.out: $(SRC_DIR)/%.c $(BIN_DIR) $(LIB_DIR)/libopenmmpol.so
	$(CC) $(CFLAGS) $(LDLIBS) -L$(LIB_DIR) -I$(INCLUDE_DIR) $< -lopenmmpol -o $@

$(LIB_DIR)/mmpolmodules.a: $(OBJS) $(LIB_DIR)
	ar crs $@ $(OBJS)

$(LIB_DIR)/libopenmmpol.so: $(OBJS) $(LIB_DIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) $(LDLIBS) -shared -I$(MOD_DIR) -o $@ $(OBJS)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)
$(MOD_DIR):
	@mkdir -p $(MOD_DIR)
$(BIN_DIR):
	@mkdir -p $(BIN_DIR)
$(LIB_DIR):
	@mkdir -p $(LIB_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f03 $(OBJ_DIR) $(MOD_DIR)
	$(FC) $(CPPFLAGS) $(FFLAGS) -I$(MOD_DIR) -J$(MOD_DIR) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR) $(MOD_DIR) $(LIB_DIR) $(DOC_DIR)

# Explicit dependencies
$(OBJ_DIR)/coulomb_kernel.o: $(OBJ_DIR)/mod_mmpol.o
$(OBJ_DIR)/elstat.o:
$(OBJ_DIR)/electrostatics.o: $(OBJ_DIR)/elstat.o $(OBJ_DIR)/mod_mmpol.o $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/energy.o: $(OBJ_DIR)/mod_mmpol.o
$(OBJ_DIR)/mmpol_init.o: $(OBJ_DIR)/mod_mmpol.o $(OBJ_DIR)/mod_memory.o 
$(OBJ_DIR)/mod_constants.o: $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/mod_interface.o: $(OBJ_DIR)/mod_mmpol.o $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/mod_io.o: $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/mod_memory.o:
$(OBJ_DIR)/mod_mmpol.o: $(OBJ_DIR)/mod_memory.o $(OBJ_DIR)/mod_constants.o $(OBJ_DIR)/mod_io.o
$(OBJ_DIR)/multipoles_functions.o: $(OBJ_DIR)/elstat.o 
$(OBJ_DIR)/polar.o: $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/polarization.o: $(OBJ_DIR)/solvers.o $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/rotate_multipoles.o: $(OBJ_DIR)/mod_mmpol.o 
$(OBJ_DIR)/solvers.o: $(OBJ_DIR)/mod_memory.o
$(OBJ_DIR)/utilities.o: $(OBJ_DIR)/mod_mmpol.o
