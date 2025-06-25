#!/bin/bash

export SUPERLU_SRC_DIR=$PWD/SuperLU_4.3
export SUPERLU_INSTALL_DIR=$HOME/local/gcc/SuperLU


cd $SUPERLU_SRC_DIR
mkdir shared
cd shared

make -f ../SRC/Makefile VPATH=../SRC srcdir=../SRC CC=cc CFLAGS="$CFLAGS -fPIC" FORTRAN=gfortran FFLAGS="$CFLAGS \
     -fPIC" PLAT="" BLASDEF="" BLASLIB="-L$HOME/local/gcc/lapack/lib64 -lblas" CDEFS="-DAdd_" NOOPTS="-fPIC" \
     ARCH="echo" ARCHFLAGS="" RANLIB="echo" SUPERLULIB=../lib/libsuperlu.a

gcc -shared -Wl,-soname,libsuperlu.so.4 -o ../lib/libsuperlu.so.4 *.o -lblas -lm -lgfortran

if [ ! -d "$SUPERLU_INSTALL_DIR/include" ]; then
   mkdir -p $SUPERLU_INSTALL_DIR/include
fi

if [ ! -d "$SUPERLU_INSTALL_DIR/lib" ]; then
   mkdir $SUPERLU_INSTALL_DIR/lib
fi


cp $SUPERLU_SRC_DIR/lib/libsuperlu.so.4 $SUPERLU_INSTALL_DIR/lib
cp $SUPERLU_SRC_DIR/SRC/*.h $SUPERLU_INSTALL_DIR/include

echo "In case of compilation problems with Trilinos, consider"
echo "commenting out the line"
echo 'extern int     lsame_ (char *, char *);'
echo "somewhere around line 349 in the header file"
echo "slu_util.h"

