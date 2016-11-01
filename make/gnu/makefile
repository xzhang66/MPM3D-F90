mpm3d-f90 :  Constitution.o  DataIn.o DataOut.o FFI.o \
         Grid.o MPM3D.o Material.o Particle.o update_step.o
	gfortran -o  mpm3d-f90 Constitution.o DataIn.o DataOut.o FFI.o \
                     Grid.o MPM3D.o  Material.o Particle.o update_step.o

Constitution.o : ../../src/Constitution.f90 MaterialData.mod ParticleData.mod FFI.mod GridData.mod 
	gfortran -c  ../../src/Constitution.f90

DataIn.o : ../../src/DataIn.f90 DataOut.mod
	gfortran -c ../../src/DataIn.f90

DataOut.mod : ../../src/DataOut.f90
	gfortran -c ../../src/DataOut.f90

DataOut.o : ../../src/DataOut.f90
	gfortran -c ../../src/DataOut.f90

FFI.mod : ../../src/FFI.F90
	gfortran -c ../../src/FFI.f90

FFI.o : ../../src/FFI.F90
	gfortran -c ../../src/FFI.f90

GridData.mod : ../../src/Grid.f90
	gfortran -c ../../src/Grid.f90

Grid.o : ../../src/Grid.f90
	gfortran -c ../../src/Grid.f90

MPM3D.o : ../../src/MPM3D.f90
	gfortran -c ../../src/MPM3D.f90

MaterialData.mod : ../../src/Material.f90
	gfortran -c ../../src/Material.f90

Material.o : ../../src/Material.f90
	gfortran -c ../../src/Material.f90

ParticleData.mod : ../../src/Particle.f90
	gfortran -c ../../src/Particle.f90

Particle.o : ../../src/Particle.f90
	gfortran -c ../../src/Particle.f90

update_step.o : ../../src/update_step.f90
	gfortran -c ../../src/update_step.f90

clean : 
	rm mpm3d-f90 *.o *.mod 
