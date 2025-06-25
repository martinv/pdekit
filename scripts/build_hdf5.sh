#!/bin/bash

# TO AVOID CRASHES, BUILD THE PARALLEL VERSION OF HDF5 __ONLY__ AS STATIC!!!

#export COMPILER_PATH=/home/martin/scfd_deps/openmpi/bin
#export CC="$COMPILER_PATH/mpicc"
#export CXX="$COMPILER_PATH/mpicxx"

#            --enable-cxx=yes \
#            --with-default-api-version=v16
#            --enable-shared=yes \
#            --enable-static=yes \

./configure --prefix=/home/martin/local/scfd_deps/hdf5-1.8 \
            --enable-parallel=yes \
            --enable-hl=yes
