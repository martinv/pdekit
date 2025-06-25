# ==============================================================================
# User-defined options for compilation
# ==============================================================================

#  set(BaseName "binary" CACHE STRING "BaseName chosen by the user at CMake configure time")

set(PDEKIT_MPI_Fortran_COMPILER /usr/bin/mpif90 CACHE FILEPATH "MPI Fortran compiler executable.")
set(PDEKIT_MPI_C_COMPILER /usr/bin/mpicc CACHE FILEPATH "MPI C compiler executable.")
set(PDEKIT_MPI_CXX_COMPILER /usr/bin/mpicxx CACHE FILEPATH "MPI C++ compiler executable.")

set(PDEKIT_BLAS_LAPACK_PROVIDER "Netlib" CACHE STRING "Select the provider of Blas and Lapack libraries. Available options: Netlib, OpenBLAS, IntelMKL, Internal (will build OpenBLAS)")
# Now, after defining the cache entry with its initial default value,
# define the set of strings to which its value should be constrained:
set_property(CACHE PDEKIT_BLAS_LAPACK_PROVIDER PROPERTY STRINGS Netlib OpenBLAS IntelMKL Internal)

set(PDEKIT_BLAS_DIR /usr/lib CACHE PATH "Directory containing Blas libraries")
set(PDEKIT_LAPACK_DIR /usr/lib CACHE PATH "Directory containing Lapack libraries")

option(PDEKIT_USE_TRILINOS "Use the Trilinos library for linear algebra" OFF)
set(PDEKIT_TRILINOS_DIR /usr/lib CACHE PATH "Directory containing the Trilinos library")
set(BOOST_ROOT /usr/lib CACHE PATH "Directory containing Boost C++ libraries")



