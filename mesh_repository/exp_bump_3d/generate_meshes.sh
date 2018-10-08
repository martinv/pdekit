#!/bin/bash

 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/opencascade/lib
 export PETSC_DIR=$HOME/local/gcc/petsc-3.5.3
 export PETSC_ARCH=linux-gnu-cxx-opt
 export SLEPC_DIR=$HOME/local/gcc/slepc-3.5.3
 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib:$SLEPC_DIR/$PETSC_ARCH/lib

 export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/gcc/parmetis/lib:$HOME/local/gcc/cgns/lib:$HOME/local/gcc/hdf5-1.8.14/lib

 export PYTHONPATH=$PYTHONPATH:/data/software/gmsh-svn/projects/dg/build/gmsh:/data/software/gmsh-svn/projects/dg/build
 export DG_BUILD_DIR=/data/software/gmsh-svn/projects/dg/build
 export MPILIBDIR=$HOME/local/gcc/openmpi/lib

 export PYTHONPATH=$PYTHONPATH:/data/software/gmsh-svn/build:/data/software/gmsh-svn/build/wrappers
 export LD_PRELOAD=:$HOME/local/gcc/openmpi/lib/libmpi_cxx.so:$HOME/local/gcc/openmpi/lib/libmpi_usempif08.so:$HOME/local/gcc/openmpi/lib/libmpi.so:$HOME/local/gcc/openmpi/lib/libopen-pal.so:$HOME/local/gcc/openmpi/lib/libopen-rte.so

gmsh -3 exp_bump_p1_tet.geo -order 1 -optimize -o exp_bump_p1_tet.msh
gmsh -3 exp_bump_p2_tet.geo -order 2 -optimize_ho -ho_min 0.7 -o exp_bump_p2_tet.msh
gmsh -3 exp_bump_p3_tet.geo -order 3 -optimize_ho -ho_min 0.7 -o exp_bump_p3_tet.msh
gmsh -3 exp_bump_p4_tet.geo -order 4 -optimize_ho -ho_min 0.7 -o exp_bump_p4_tet.msh

