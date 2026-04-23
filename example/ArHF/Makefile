# Makefile for IALR_WP_CRWP (scatt)

# Compiler and flags
FC = gfortran

# Base flags for production
FFLAGS_BASE = -O3 -fopenmp

# Optional extra speed when FAST=1 (CPU-specific + aggressive math)
FFASTFLAGS = -march=native -funroll-loops -ffast-math

# Effective flags
FFLAGS = $(FFLAGS_BASE)
ifeq ($(FAST),1)
FFLAGS += $(FFASTFLAGS)
endif

# Debug-only extras; keep off for production runs because they add large overhead
FFLAGS_STRICT = $(FFLAGS_BASE) -fcheck=bounds -fbacktrace -g 
#FFLAGS_STRICT = $(FFLAGS_BASE) -fcheck=bounds -fbacktrace -g -fcheck=bounds
FFLAGS_LEGACY = $(FFLAGS) -std=legacy -fallow-argument-mismatch
LDFLAGS = -fopenmp -framework Accelerate

# Target executable
TARGET = scatt
CS_TARGET = CrossSection

# Source files
# Fixed format Fortran files
FSRC = dfftpack.f function.f interaction-PES.f ArHF_PES.f

# Free format Fortran 90 files (in compilation order based on dependencies)
F90SRC = machina_basic.f90 \
        globalPara.f90 \
        potent.f90 \
        basis.f90 \
        initWP.f90 \
        interCoorRCB.f90 \
        propagation.f90 \
        match.f90 \
        main.f90

# Object files
FOBJ = $(FSRC:.f=.o)
F90OBJ = $(F90SRC:.f90=.o)
OBJ = $(FOBJ) $(F90OBJ)
CS_OBJ = CrossSection.o

# Default target
all: $(TARGET) $(CS_TARGET)


# Link object files to create executable
$(TARGET): $(OBJ)
	$(FC) $(LDFLAGS) -o $@ $(OBJ)
	@echo "Successfully built $(TARGET)"

# Build CrossSection executable (needs function.f and m_MachinaBasic from machina_basic.f90)
$(CS_TARGET): machina_basic.o function.o $(CS_OBJ)
	$(FC) $(LDFLAGS) -o $@ $^
	@echo "Successfully built $(CS_TARGET)"

# Compile fixed format Fortran files (legacy code, allow old-style constructs)
%.o: %.f
	$(FC) $(FFLAGS_LEGACY) -c $<

# Compile free format Fortran 90 files (fast by default; use `make DEBUG=1` for checks)
%.o: %.f90
ifeq ($(DEBUG),1)
	$(FC) $(FFLAGS_STRICT) -c $<
else
	$(FC) $(FFLAGS) -c $<
endif

# Module dependencies (ensure modules are compiled before files that use them)
globalPara.o: machina_basic.o
potent.o: machina_basic.o globalPara.o
basis.o: machina_basic.o globalPara.o potent.o
initWP.o: machina_basic.o globalPara.o basis.o
interCoorRCB.o: machina_basic.o globalPara.o basis.o potent.o
propagation.o: machina_basic.o globalPara.o initWP.o potent.o basis.o
match.o: machina_basic.o globalPara.o interCoorRCB.o
main.o: machina_basic.o globalPara.o potent.o basis.o initWP.o propagation.o interCoorRCB.o match.o
CrossSection.o: machina_basic.o

# Clean up
clean:
	rm -f $(OBJ) *.mod

distclean: clean
	rm -f $(TARGET) $(CS_TARGET)

# Phony targets
.PHONY: all clean distclean
