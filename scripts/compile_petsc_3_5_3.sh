#!/bin/bash

export PETSC_DIR=$PWD
export PETSC_ARCH="linux-gnu-cxx-opt"

# ==============================================================================

# TO COMPILE PETSC 3.5.3
python2 ./config/configure.py --COPTFLAGS=O2 --CXXOPTFLAGS=-O2 --CFLAGS="-O2" --CXXFLAGS="-O2 -fpermissive" -with-clanguage=cxx --with-c++-support --with-shared-libraries=1 --with-numpy --with-blas-lapack-lib-dir=$HOME/local/lapack/lib --download-metis=1 --download-parmetis=1 --download-suitesparse=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 --download-blacs=1 
#--download-pastix=1
#--download-ptscotch=1 

# ==============================================================================

make 

