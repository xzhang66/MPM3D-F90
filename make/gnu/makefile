# Compiler flags
FFLAGS ?= -O2
FFLAGS_DEBUG = -g -O0 -fcheck=all -fbacktrace -Wall -Wextra
FFLAGS_RELEASE = -O3 -ffast-math -march=native -DNDEBUG

# Default target
all: mpm3d-f90

# Release mode
release: FFLAGS = $(FFLAGS_RELEASE)
release: clean mpm3d-f90

# Debug mode  
debug: FFLAGS = $(FFLAGS_DEBUG)
debug: clean mpm3d-f90

# Main executable
mpm3d-f90 : Constitution.o DataIn.o DataOut.o FFI.o \
         Grid.o MPM3D.o Material.o Particle.o update_step.o
	gfortran $(FFLAGS) -o mpm3d-f90 Constitution.o DataIn.o DataOut.o FFI.o \
                     Grid.o MPM3D.o Material.o Particle.o update_step.o

# Object file rules with FFLAGS
Constitution.o : ../../src/Constitution.f90 MaterialData.mod ParticleData.mod FFI.mod GridData.mod
	gfortran $(FFLAGS) -c ../../src/Constitution.f90

DataIn.o : ../../src/DataIn.f90 DataOut.mod
	gfortran $(FFLAGS) -c ../../src/DataIn.f90

DataOut.mod : ../../src/DataOut.f90
	gfortran $(FFLAGS) -c ../../src/DataOut.f90

DataOut.o : ../../src/DataOut.f90
	gfortran $(FFLAGS) -c ../../src/DataOut.f90

FFI.mod : ../../src/FFI.f90
	gfortran $(FFLAGS) -c ../../src/FFI.f90

FFI.o : ../../src/FFI.f90
	gfortran $(FFLAGS) -c ../../src/FFI.f90

GridData.mod : ../../src/Grid.f90
	gfortran $(FFLAGS) -c ../../src/Grid.f90

Grid.o : ../../src/Grid.f90
	gfortran $(FFLAGS) -c ../../src/Grid.f90

MPM3D.o : ../../src/MPM3D.f90
	gfortran $(FFLAGS) -c ../../src/MPM3D.f90

MaterialData.mod : ../../src/Material.f90
	gfortran $(FFLAGS) -c ../../src/Material.f90

Material.o : ../../src/Material.f90
	gfortran $(FFLAGS) -c ../../src/Material.f90

ParticleData.mod : ../../src/Particle.f90
	gfortran $(FFLAGS) -c ../../src/Particle.f90

Particle.o : ../../src/Particle.f90
	gfortran $(FFLAGS) -c ../../src/Particle.f90

update_step.o : ../../src/update_step.f90
	gfortran $(FFLAGS) -c ../../src/update_step.f90

# Utility targets
clean:
	rm -f mpm3d-f90 *.o *.mod

help:
	@echo "Available targets:"
	@echo "  make all       - Build with default optimization (-O2)"
	@echo "  make release   - Build with maximum optimization (-O3)"
	@echo "  make debug     - Build with debug flags"
	@echo "  make clean     - Remove all compiled files"

.PHONY: all release debug clean help