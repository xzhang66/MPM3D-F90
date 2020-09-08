cmake_minimum_required(VERSION 3.2)

project(mpm3d)
enable_language(Fortran)

file(GLOB_RECURSE sources  src/*.f90 src/*.h)
add_executable(mpm3d ${sources})
