# ====================================================================

EXTRA_ARGS=$@

# ====================================================================
# USER-DEFINED VARIABLES
#
export MY_TRILINOS_DEPS_DIR=$HOME/local/gcc
export MY_TRILINOS_INSTALL_DIR=$HOME/local/gcc/trilinos-12.6.3
export MY_TRILINOS_SRC_DIR=$PWD/trilinos-12.6.3-Source
export MY_TRILINOS_BUILD_DIR=$MY_TRILINOS_SRC_DIR/build
#
# ====================================================================

cd $MY_TRILINOS_BUILD_DIR

# this causes build failures on my workstation
#  -D CMAKE_BUILD_TYPE:STRING=DEBUG \

# To include HDF5 support:
#  -D EpetraExt_ENABLE_HDF5:BOOL=ON \
#  -D TPL_ENABLE_XDMF:BOOL=ON \

#  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON \

cmake \
  -D CMAKE_Fortran_COMPILER=$MY_TRILINOS_DEPS_DIR/openmpi/bin/mpif90 \
  -D CMAKE_C_COMPILER=$MY_TRILINOS_DEPS_DIR/openmpi/bin/mpicc \
  -D CMAKE_CXX_COMPILER=$MY_TRILINOS_DEPS_DIR/openmpi/bin/mpicxx \
  -D CMAKE_EXE_LINKER_FLAGS="-llapack -lblas" \
  -D MPI_EXEC=$MY_TRILINOS_DEPS_DIR/openmpi/bin/mpiexec \ -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D BUILD_SHARED_LIBS:BOOL=ON \
  -D Trilinos_ENABLE_CXX11=ON \
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D ML_ENABLE_Aztec:BOOL=ON \
  -D ML_ENABLE_EXAMPLES:BOOL=ON \
  -D ML_ENABLE_Enrich:BOOL=ON \
  -D ML_ENABLE_MLapi:BOOL=ON \
  -D ML_ENABLE_MPI:BOOL=ON \
  -D ML_ENABLE_SUPERLUDIST:BOOL=ON \
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Trilinos_ENABLE_Amesos2:BOOL=ON \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_Belos:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
  -D Trilinos_ENABLE_MueLu:BOOL=ON \
  -D MueLu_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_NOX:BOOL=ON \
  -D Trilinos_ENABLE_Teko:BOOL=ON \
  -D Trilinos_ENABLE_Rythmos:BOOL=ON \
  -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
  -D Trilinos_ENABLE_Sundance:BOOL=OFF \
  -D Trilinos_ENABLE_Thyra:BOOL=ON \
  -D Trilinos_ENABLE_Tpetra:BOOL=ON \
  -D Trilinos_ENABLE_Xpetra:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
  -D Zoltan2_ENABLE_EXAMPLES:BOOL=ON \
  -D Anasazi_ENABLE_TESTS:BOOL=ON \
  -D Anasazi_ENABLE_EXAMPLES:BOOL=ON \
  -D Amesos2_ENABLE_Basker:BOOL=ON \
  -D Amesos2_ENABLE_KLU2:BOOL=ON \
  -D EpetraExt_ENABLE_AMD:BOOL=ON \
  -D Epetra_ENABLE_TESTS:BOOL=OFF \
  -D Epetra_ENABLE_EXAMPLES:BOOL=OFF \
  -D DART_TESTING_TIMEOUT:STRING=600 \
  -D TPL_ENABLE_MPI:BOOL=ON \
  -D TPL_ENABLE_HDF5:BOOL=ON \
  -D HDF5_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/hdf5-1.8.16/include \
  -D HDF5_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/hdf5-1.8.16/lib \
  -D AMD_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/suitesparse/include \
  -D AMD_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/suitesparse/lib \
  -D TPL_ENABLE_BLAS:BOOL=ON \
  -D BLAS_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/lapack/lib64 \
  -D BLAS_LIBRARY_NAMES:STRING="blas" \
  -D TPL_ENABLE_LAPACK:BOOL=ON \
  -D LAPACK_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/lapack/lib64 \
  -D TPL_ENABLE_UMFPACK:BOOL=ON \
  -D UMFPACK_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/suitesparse/include \
  -D UMFPACK_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/suitesparse/lib \
  -D UMFPACK_LIBRARY_NAMES:STRING="umfpack;amd;cholmod;camd;colamd;ccolamd;ufconfig" \
  -D TPL_ENABLE_METIS:BOOL=ON \
  -D METIS_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/include \
  -D METIS_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/lib \
  -D TPL_ENABLE_ParMETIS:BOOL=ON \
  -D ParMETIS_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/include \
  -D ParMETIS_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/lib \
  -D TPL_ENABLE_Boost:BOOL=ON \
  -D Boost_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/boost/include \
  -D Boost_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/boost/lib \
  -D TPL_ENABLE_Netcdf:BOOL=OFF \
  -D TPL_ENABLE_SuperLU:BOOL=ON \
  -D SuperLU_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/SuperLU/include \
  -D SuperLU_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/SuperLU/lib \
  -D TPL_SuperLU_LIBRARIES:PATH="${MY_TRILINOS_DEPS_DIR}/SuperLU/lib/libsuperlu.so.4" \
  -D TPL_ENABLE_Pthread:BOOL=ON \
  -D TPL_ENABLE_Zlib:BOOL=ON \
  -D CMAKE_INSTALL_PREFIX:PATH=${MY_TRILINOS_INSTALL_DIR} \
  $EXTRA_ARGS \
  ${MY_TRILINOS_SRC_DIR}



#  -D METIS_INCLUDE_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/include \
#  -D METIS_LIBRARY_DIRS:PATH=${MY_TRILINOS_DEPS_DIR}/parmetis/lib \
#  -D Zoltan2_ENABLE_METIS:BOOL=OFF \
#  -D Zoltan2_ENABLE_ParMETIS:BOOL=OFF \
#  -D Zoltan2_ENABLE_ANASAZI:BOOL=OFF \








  #-D UMFPACK_LIBRARY_NAMES:STRING="umfpack;amd;cholmod;camd;colamd;ccolamd;ufconfig;metis" \
  #-D EpetraExt_ENABLE_Zoltan:BOOL=ON \
  #-D EpetraExt_BUILD_Zoltan:BOOL=ON \
  #-D EpetraExt_BUILD_PETSC:BOOL=ON \
  #-D EpetraExt_BUILD_UMFPACK:BOOL=ON \
  #-D ML_ENABLE_Zoltan2:BOOL=ON \
  #-D Trilinos_Anasazi_ENABLE_Epetra:BOOL=ON \
