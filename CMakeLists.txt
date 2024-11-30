cmake_minimum_required(VERSION 3.16)

project(Paraprog)

set(CXX_STANDARD 17)
set(CXX_STANDARD_REQUIRED On)

option (ADD_MPI "Add MPI library for parallel matrix library" ON)

#Find MPI
if(ADD_MPI)
    find_package(MPI)
endif()

if(MPI_FOUND)
    set(MPI_C_COMPILER "/usr/share/mvapich/bin/mpicc")
    add_compile_definitions("USE_MPI")
endif()

add_subdirectory(MatrixLib)
add_subdirectory(Main)