#!/usr/bin/python2

import time
import sys
from subprocess import call
from optparse import OptionParser
import os

# ========================================================================================================
# Python script for installation of various software dependencies. See install_scfd_deps.py -h for usage
# ========================================================================================================
# Compiler variables: CC, CXX, F77, or F90
#
#
#
#
# ========================================================================================================

def install_cmake(optlist):
  print "Starting installation of cmake ... "
  
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd="tar -xzf cmake-3.2.1.tar.gz"
  call(cmnd,shell=True)
  
  cmake_src_dir     = os.path.abspath(optlist.repodir + "/cmake-3.2.1")
  
  if ( os.path.exists(optlist.cmake_install_dir) ):
    print "The folder %s already exists. Will not install. Exiting." % cmake_install_dir
    exit(1)
    
  os.chdir(cmake_src_dir)
  
  cmnd = "./bootstrap --prefix=%s" % optlist.cmake_install_dir
  call(cmnd,shell=True)
  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "make install"
  call(cmnd,shell=True)
  
  os.chdir(rootpath)
  
  print " ... cmake installed"


# ========================================================================================================

def install_openmpi(optlist):
  print "Starting installation of openmpi ... "
  
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd="tar -xzf openmpi-3.1.0.tar.gz"
  call(cmnd,shell=True)
  
  openmpi_src_dir     = os.path.abspath(optlist.repodir + "/openmpi-3.1.0")

  if ( os.path.exists(optlist.mpi_install_dir) ):
    print "The folder %s already exists. Will not install openmpi. Exiting." % openmpi_install_dir
    exit(1)
  
  os.chdir(openmpi_src_dir)

  # CC=icc CXX=icpc F77=ifort FC=ifort 
  #cmnd = "./configure --enable-shared --enable-static --enable-mpi-threads --with-threads=posix \
  #--prefix=" + openmpi_install_dir + " --with-mpi-f90-size=medium"
  #cmnd = "./configure CC=%s CXX=%s F77=%s FC=%s --enable-shared --enable-mpi-threads --with-threads=posix \
  #--prefix=" + openmpi_install_dir + " --with-mpi-f90-size=medium" % (optlist.ccompiler,optlist.cppcompiler,\
  #                                                                    optlist.fortrancompiler,optlist.fortrancompiler)
  cmnd = "./configure CC=" + optlist.ccompiler + " CXX=" + optlist.cppcompiler +\
         " FC=" + optlist.fortrancompiler +\
         " --enable-shared --enable-mpi-thread-multiple --with-threads=posix \
           --prefix=" + optlist.mpi_install_dir + " --with-mpi-f90-size=medium"
  #cmnd = "./configure CC=" + optlist.ccompiler + " CXX=" + optlist.cppcompiler +\
  #       " F77=" + optlist.fortrancompiler + " FC=" + optlist.fortrancompiler +\
  #       " --prefix=" + optlist.mpi_install_dir

  call(cmnd,shell=True)
  call(optlist.makecmnd,shell=True)
  cmnd="make install"
  call(cmnd,shell=True)
  
  os.chdir(rootpath)
  
  print " ... openmpi installed"
  
# ========================================================================================================

def install_lapack(optlist):
  # Note: lapack svn repository is in 
  # svn co https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk
  print "Starting installation of lapack ... "
  
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd="tar -xzf lapack-3.8.0.tar.gz"
  call(cmnd,shell=True)
  
  lapack_src_dir     = os.path.abspath(optlist.repodir + "/lapack-3.8.0")
  lapack_build_dir   = os.path.abspath(optlist.repodir + "/lapack-3.8.0/build")

  if ( os.path.exists(optlist.lapack_install_dir) ):
    print "The folder %s already exists. Will not install lapack. Exiting." % optlist.lapack_install_dir
    exit(1)
 
  os.mkdir( lapack_build_dir )
 
  os.chdir(lapack_build_dir)
  
 #-D BUILD_STATIC_LIBS:BOOL=ON \

  print "FORTRAN COMPILER: %s" % optlist.fortrancompiler

  cmnd = "cmake -D CMAKE_Fortran_COMPILER=" + optlist.fortrancompiler + " \
 -D BUILD_SHARED_LIBS:BOOL=ON \
 -D CMAKE_BUILD_TYPE:STRING=Release \
 -D USE_OPTIMIZED_BLAS:BOOL=OFF \
 -D LAPACKE:BOOL=ON \
 -D CMAKE_INSTALL_PREFIX:PATH=" + optlist.lapack_install_dir  + " " + lapack_src_dir

 # -D LAPACKE_WITH_TMG:BOOL=ON \
  
  call(cmnd,shell=True)

  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "make install"
  call(cmnd,shell=True)

  os.chdir(rootpath)

  cmnd="rm -rf lapack-3.8.0"
  call(cmnd,shell=True)
  
  print " ... lapack installed"


# ========================================================================================================

def install_atlas(optlist):
  
  print "Starting installation of atlas ..."
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  atlas_src_dir     = os.path.abspath(optlist.repodir+"/ATLAS")
  atlas_build_dir   = os.path.abspath(optlist.repodir+"/ATLAS/build")
  atlas_install_dir = os.path.abspath(optlist.installdir + "/atlas")

  if ( os.path.exists(atlas_install_dir) ):
    print "The folder %s already exists. Will not install atlas. Exiting." % atlas_install_dir
    exit(1)

  cmnd = "tar -xjf atlas3.9.45.tar.bz2"
  call(cmnd,shell=True)
  os.mkdir(atlas_build_dir)
  os.chdir(atlas_build_dir)

  lapack_archive = os.path.abspath(optlist.repodir + "/lapack-3.3.1.tgz")

  #cmnd = "../configure --prefix=%s -shared -Fa alg -fPIC -Si cputhrchk 0 --with-netlib-lapack-tarfile=%s" % (atlas_install_dir,lapack_archive)
  cmnd = "../configure --prefix=%s --shared -Fa alg -fPIC --with-netlib-lapack-tarfile=%s" % (atlas_install_dir,lapack_archive)
  print "The command is "
  print cmnd
  call(cmnd,shell=True)
  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "make install"
  call(cmnd,shell=True)

  os.chdir(rootpath)
  
  print " ... atlas installed"
  

# ========================================================================================================

def install_suitesparse(optlist):
  print "Starting installation of suitesparse (amd,umfpack) ..."

  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  #cmnd = "tar -xzf SuiteSparse-3.6.1-TintinMod.tar.gz"
  cmnd = "tar -xzf SuiteSparse.tar.gz"
  call(cmnd,shell=True)

  #suitesparse_src_dir   = os.path.abspath(optlist.repodir+"/SuiteSparse-TintinMod")
  suitesparse_src_dir   = os.path.abspath(optlist.repodir+"/SuiteSparse")

  suitesparse_install_dir         = os.path.abspath(optlist.installdir + "/suitesparse")
  suitesparse_install_include_dir = os.path.abspath(optlist.installdir + "/suitesparse/include")
  suitesparse_install_lib_dir      = os.path.abspath(optlist.installdir + "/suitesparse/lib")

  if ( os.path.exists(suitesparse_install_dir) ):
    print "The folder %s already exists. Will not install suitesparse. Exiting." % suitesparse_install_dir
    exit(1)


  ufconfig_dir = os.path.abspath(suitesparse_src_dir + "/UFconfig") 
  os.chdir(ufconfig_dir)

  cmnd = "mv UFconfig.mk backup_UFconfig.mk"
  call(cmnd,shell=True)

  # Create new customized file UFconfig.mk
  outfile = open("UFconfig.mk",'w')

  outputstr = """#===============================================================================
# UFconfig.mk:  common configuration file for the SuiteSparse
#===============================================================================

#------------------------------------------------------------------------------
# Generic configuration
#------------------------------------------------------------------------------

# Using standard definitions from the make environment, typically:
#
#   CC              cc      C compiler
#   CXX             g++     C++ compiler
#   CFLAGS          [ ]     flags for C and C++ compiler
#   CPPFLAGS        [ ]     flags for C and C++ compiler
#   TARGET_ARCH     [ ]     target architecture
#   FFLAGS          [ ]     flags for Fortran compiler
#   RM              rm -f   delete a file
#   AR              ar      create a static *.a library archive
#   ARFLAGS         rv      flags for ar
#   MAKE            make    make itself (sometimes called gmake)
#
# You can redefine them here, but by default they are used from the
# default make environment.

# C and C++ compiler flags.  The first three are standard for *.c and *.cpp
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O3 -fexceptions -fPIC

# ranlib, and ar, for generating libraries. If you don't need ranlib,
# just change it to RANLAB = echo
RANLIB = ranlib
ARCHIVE = $(AR) $(ARFLAGS)

# copy, delete, and rename a file
CP = cp -f
MV = mv -f
# RM = rm -f

# Fortran compiler (not required for 'make' or 'make library')
F77 = gfortran
F77FLAGS = $(FFLAGS) -O
F77LIB =

# C and Fortran libraries
LIB = -lm

# For compiling MATLAB mexFunctions (MATLAB 7.5 or later)
MEX = mex -O -largeArrayDims -lmwlapack -lmwblas

# For compiling MATLAB mexFunctions (MATLAB 7.3 and 7.4)
# MEX = mex -O -largeArrayDims -lmwlapack

# For MATLAB 7.2 or earlier, you must use one of these options:
# MEX = mex -O -lmwlapack
# MEX = mex -O

# Which version of MAKE you are using (default is "make")
# MAKE = make
# MAKE = gmake

# For "make install"
"""

  outfile.write(outputstr)

  outputstr = "INSTALL_LIB = %s\n" % suitesparse_install_lib_dir
  outfile.write(outputstr)
  outputstr = "INSTALL_INCLUDE = %s\n\n" % suitesparse_install_include_dir
  outfile.write(outputstr)

  outputstr = """#------------------------------------------------------------------------------
# BLAS and LAPACK configuration:
#------------------------------------------------------------------------------
"""

  outfile.write(outputstr)
  
  #outputstr = "BLAS   = -L%s -llapack -lf77blas -lcblas -latlas -lgfortran\n" % os.path.abspath(optlist.depsroot + "/atlas/lib")
  #outputstr = "BLAS   = -L%s -lblas -llapack -lgfortran\n" % os.path.abspath(optlist.depsroot + "/lapack/lib")
  #outfile.write(outputstr)
  outputstr = "BLAS   = -L%s -lblas -llapack -lgfortran\n" % os.path.abspath(optlist.lapack_lib_dir)
  outfile.write(outputstr)  
  #outputstr = "LAPACK = -L%s -llapack -lf77blas -lcblas -latlas -lgfortran\n" % os.path.abspath(optlist.depsroot + "/atlas/lib")
  #outputstr = "LAPACK = -L%s -lblas -llapack -lgfortran\n" % os.path.abspath(optlist.depsroot + "/lapack/lib")
  outputstr = "LAPACK = -L%s -lblas -llapack -lgfortran\n" % os.path.abspath(optlist.lapack_lib_dir)
  outfile.write(outputstr)
  outputstr = "XERBLA = \n\n"
  outfile.write(outputstr)


  outputstr = """#------------------------------------------------------------------------------
# METIS, optionally used by CHOLMOD
#------------------------------------------------------------------------------
"""
  outfile.write(outputstr)

#METIS_PATH = "../../metis-4.0.3
#METIS = ../../metis-4.0.3/libmetis.a

  outputstr = "METIS_PATH = " + optlist.parmetis_dir + "\nMETIS = " + optlist.parmetis_lib_dir + "/libmetis.a\n"
  outfile.write(outputstr)

  outputstr="""
#------------------------------------------------------------------------------
# UMFPACK configuration:
#------------------------------------------------------------------------------

# Configuration flags for UMFPACK.  See UMFPACK/Source/umf_config.h for details.

UMFPACK_CONFIG =

#------------------------------------------------------------------------------
# CHOLMOD configuration
#------------------------------------------------------------------------------


CHOLMOD_CONFIG =

#------------------------------------------------------------------------------
# SuiteSparseQR configuration:
#------------------------------------------------------------------------------

# default, without timing, without TBB:
SPQR_CONFIG =

# with TBB, you must select this:
# TBB = -ltbb
# without TBB:
TBB =

# with timing, you must include the timing library:
# RTLIB = -lrt
# without timing
RTLIB =

#------------------------------------------------------------------------------
# remove object files and profile output
#------------------------------------------------------------------------------

CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d *.gcda *.gcno

"""
  outfile.write(outputstr)
  outfile.close()

  os.mkdir(suitesparse_install_dir)
  os.mkdir(suitesparse_install_include_dir)
  os.mkdir(suitesparse_install_lib_dir)

  # Compile UFconfig:
  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "make install"
  call(cmnd,shell=True)

  """
  # Compile and install metis 4.0.3:
  metisdir = os.path.abspath(suitesparse_src_dir + "/metis-4.0.3")
  os.chdir(metisdir) 
  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "cp Lib/metis.h %s" % suitesparse_install_include_dir
  call(cmnd,shell=True)
  cmnd = "cp libmetis.a %s" % suitesparse_install_lib_dir
  call(cmnd,shell=True)

  # Compile all packages needed by umfpack and umfpack itself. Install all of them
  for package in ["AMD","COLAMD","CCOLAMD","CAMD","CHOLMOD","UMFPACK"]:
    workdir = os.path.abspath(suitesparse_src_dir + "/" + package)
    print "\t Suitesparse: compiling package %s " % package
    os.chdir(workdir)
    cmnd = "make"
    call(cmnd,shell=True)
    cmnd = "make install"
    call(cmnd,shell=True)
  """


  os.chdir(suitesparse_src_dir)
  cmnd = "make"
  call(cmnd,shell=True)
  cmnd = "make install"
  call(cmnd,shell=True)

  os.chdir(rootpath)

  print " ... suitesparse installed"



# ========================================================================================================

def install_ptscotch(optlist):
  print "Starting installation of PT scotch ..."

  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd = "tar -xzf scotch_5.1.12b.tar.gz"
  call(cmnd,shell=True)
  
  scotch_src_dir      = os.path.abspath(optlist.repodir+"/scotch_5.1.12/src")
  scotch_makeinc_dir  = os.path.abspath(optlist.repodir+"/scotch_5.1.12/src/Make.inc")
  os.chdir(scotch_makeinc_dir)
  
  outfile = open("Makefile_scfd.inc",'w')
  
  outputstr = """EXE             =
LIB             = .a
OBJ             = .o

MAKE            = make
AR              = ar
ARFLAGS         = -ruv
CAT             = cat
CCS             = """ + optlist.ccompiler + "\n"
  outfile.write(outputstr)
  outfile.write("\n")
  #outputstr = "CCP             = %s/openmpi/bin/mpicc\n" % optlist.installdir
  outputstr = "CCP             = %s\n" % optlist.mpiccompiler
  outfile.write(outputstr)
  
  #outputstr = "CCD             = gcc -I%s/openmpi/include\n" % optlist.installdir
  outputstr = "CCD             = %s -I%s\n" % (optlist.ccompiler,optlist.mpi_include_dir)
  outfile.write(outputstr)
 
  outputstr = """CFLAGS          = -O3 -fPIC -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -DSCOTCH_PTHREAD -Drestrict=__restrict -DIDXSIZE64
CLIBFLAGS       =
LDFLAGS         = -lz -lm -lrt -lpthread
CP              = cp
LEX             = flex -Pscotchyy -olex.yy.c
LN              = ln
MKDIR           = mkdir
MV              = mv
RANLIB          = ranlib
YACC            = bison -pscotchyy -y -b y"""

  outfile.write(outputstr)
  outfile.close()
  
  os.chdir(scotch_src_dir)
  
  cmnd = "ln -s Make.inc/Makefile_scfd.inc Makefile.inc"
  call(cmnd,shell=True)
  #Build serial scotch
  cmnd = "make"
  call(cmnd,shell=True)
  #Build parallel scotch
  cmnd = "make ptscotch"
  call(cmnd,shell=True)
  
  #Scotch doesn't create the installation directory, you have to do it yourself
  os.mkdir(optlist.scotch_install_dir)
  cmnd = "make prefix=" +  optlist.scotch_install_dir + " install"
  call(cmnd,shell=True)

  os.chdir(rootpath)
  
  print " ... PT scotch installed"

# ========================================================================================================

def install_boost(optlist):
  print "Starting installation of boost ... "
  
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd="tar -xzf boost_1_63_0.tar.gz"
  call(cmnd,shell=True)
  
  boost_src_dir     = os.path.abspath(optlist.repodir + "/boost_1_63_0")
  boost_install_dir = os.path.abspath(optlist.installdir + "/boost")
  os.chdir(boost_src_dir)
 
  if optlist.cppcompiler == "icpc":
    cmnd = "./bootstrap.sh --prefix=" + boost_install_dir + " --with-toolset=intel"
    print "Preparing to build boost with toolset: intel\n"
  else:
    cmnd = "./bootstrap.sh --prefix=" + boost_install_dir + " --with-toolset=gcc"
    print "Preparing to build boost with toolset: gcc\n"

  call(cmnd,shell=True)
  cmnd="./b2 --build-dir=build toolset=gcc threading=multi --layout=system"
  #./bjam --build-dir=build toolset=gcc threading=multi -j6
  #   or ./bjam -j6 --build-dir=build toolset=gcc --build-type=complete --layout=versioned
  #   ( another option is to use threading=single)

  call(cmnd,shell=True)
  cmnd="./b2 install"
  call(cmnd,shell=True)
  
  os.chdir(rootpath)
  
  print " ... boost installed"
  
# ========================================================================================================

def install_hdf5(optlist):
  
  print "Starting installation of hdf5 ..."

  rootpath = os.getcwd()
  os.chdir(optlist.repodir)

  #mpipath = os.path.abspath(optlist.depsroot + "/local/openmpi/bin/")
  #mpicc  = mpipath + "/mpicc"
  #mpicxx = mpipath + "/mpicxx"

  cmnd  = "export CC=" + optlist.mpiccompiler
  call(cmnd,shell=True)
  cmnd  = "export CXX=" + optlist.mpicppcompiler #mpicxx
  call(cmnd,shell=True)

  cmnd="tar -xzf hdf5-1.8.14.tar.gz"
  call(cmnd,shell=True)
  
  hdf5_src_dir     = os.path.abspath(optlist.repodir + "/hdf5-1.8.14")
  hdf5_build_dir   = os.path.abspath(optlist.repodir + "/hdf5-1.8.14/build")
  #hdf5_install_dir = os.path.abspath(optlist.installdir + "/hdf5-1.8")


  os.chdir(hdf5_src_dir)

  # Needed for cmake-based installation:
  """
  os.mkdir(hdf5_build_dir)
  os.chdir(hdf5_build_dir)
  """

  if ( os.path.exists(optlist.hdf5_install_dir) ):
    print "The folder %s already exists. Will not install. Exiting." % optlist.hdf5_install_dir
    exit(1)

  # Needed for cmake-based installation:
  """
  if optlist.hdf5_use_old_api == "yes":
    cmnd = "cmake -D BUILD_SHARED_LIBS:BOOL=ON \
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D HDF5_BUILD_HL_LIB:BOOL=ON \
    -D HDF5_ENABLE_PARALLEL:BOOL=ON \
    -D HDF5_BUILD_TOOLS:BOOL=ON \
    -D HDF5_USE_16_API_DEFAULT:BOOL=ON \
    -D CMAKE_INSTALL_PREFIX:PATH=" + optlist.hdf5_install_dir + "/hdf5-1.8.14"
  else:
    cmnd = "cmake -D BUILD_SHARED_LIBS:BOOL=ON \
    -D CMAKE_C_COMPILER=" + optlist.mpiccompiler + "\
    -D CMAKE_CXX_COMPILER=" + optlist.mpicppcompiler + "\
    -D CMAKE_BUILD_TYPE:STRING=RELEASE \
    -D HDF5_BUILD_HL_LIB:BOOL=ON \
    -D HDF5_ENABLE_PARALLEL:BOOL=ON \
    -D HDF5_BUILD_TOOLS:BOOL=ON \
    -D HDF5_USE_16_API_DEFAULT:BOOL=OFF \
    -D CMAKE_INSTALL_PREFIX:PATH=" + optlist.hdf5_install_dir + "/hdf5-1.8.14"
  """


  # The --enable-shared option might be dangerous!
  """
  --enable-cxx=yes
  --with-defaul-api-version=v16
  --enable-shared=yes
  --enable-static=yes
  """

  if optlist.hdf5_use_old_api == "yes":
    cmnd = "./configure --enable-shared=yes --enable-static=yes --enable-parallel=yes --enable-hl=yes --with-default-api-version=v16 --prefix=" + optlist.hdf5_install_dir
  else:
    cmnd = "./configure --enable-shared=yes --enable-static=yes --enable-parallel=yes --enable-hl=yes --prefix=" + optlist.hdf5_install_dir
  #print "The command = %s" % cmnd

  call(cmnd,shell=True)

  cmnd = "make"
  call(cmnd,shell=True)

  #os.mkdir(hdf5_install_dir)
  
  cmnd = "make install"
  call(cmnd,shell=True)

  os.chdir(rootpath)
  
  print " ... hdf5 installed"

# ========================================================================================================  

def install_parmetis(optlist):
  
  print "Starting installation of parmetis ... "
  
  rootpath = os.getcwd()
  os.chdir(optlist.repodir)
  
  cmnd = "tar -xzf parmetis-4.0.3.tar.gz"
  call(cmnd,shell=True)
  
  parmetis_src_dir     = os.path.abspath(optlist.repodir + "/parmetis-4.0.3")
  metis_build_dir      = os.path.abspath(optlist.repodir + "/parmetis-4.0.3/metis/build")
  parmetis_build_dir   = os.path.abspath(optlist.repodir + "/parmetis-4.0.3/build")
  
  if ( os.path.exists(optlist.parmetis_install_dir) ):
    print "The folder %s already exists. Will not install. Exiting." % parmetis_install_dir
    exit(1)

    
  gklibpath="%s/metis/GKlib" % parmetis_src_dir.rstrip('/')  
    
  # Compile and install metis
  os.chdir(metis_build_dir)
  
  cmnd = "cmake -DCMAKE_C_COMPILER=" + optlist.ccompiler + " -DCMAKE_CXX_COMPILER=" + optlist.cppcompiler +\
  " -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX:PATH=" + optlist.parmetis_install_dir + \
  " -DGKLIB_PATH:PATH=" + os.path.abspath(parmetis_src_dir + "/metis/GKlib")  + " .."
  call(cmnd,shell=True)
  
  cmnd="make"
  call(cmnd,shell=True)
  cmnd="make install"
  call(cmnd,shell=True)  
  
  # Compile and install PARmetis
  os.chdir(parmetis_build_dir)
  
  cmnd = "cmake -DCMAKE_C_COMPILER=" + optlist.mpiccompiler + " -DCMAKE_CXX_COMPILER=" + optlist.mpicppcompiler + \
  " -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_INSTALL_PREFIX:PATH=" + optlist.parmetis_install_dir + \
  " -DGKLIB_PATH:PATH=" + os.path.abspath(parmetis_src_dir + "/metis/GKlib")  + \
  " -DMETIS_PATH:PATH=" + os.path.abspath(parmetis_src_dir + "/metis/") + " .."
  call(cmnd,shell=True)
  
  cmnd="make"
  call(cmnd,shell=True)
  cmnd="make install"
  call(cmnd,shell=True)
  
  os.chdir(rootpath)
  
  print " ... parmetis installed"
  
# ========================================================================================================

def install_mmg3d(optlist):
  
  print "Starting installation of mmg3d ..."

  rootpath = os.getcwd()
  os.chdir(optlist.repodir)

  cmnd="tar -xzf mmg3d4_TintinMod.tar.gz"
  call(cmnd,shell=True)
  
  mmg3d_src_dir     = os.path.abspath(optlist.repodir + "/mmg3d4_TintinMod")
  mmg3d_build_dir   = os.path.abspath(optlist.repodir + "/mmg3d4_TintinMod/build") 
  mmg3d_install_dir = os.path.abspath(optlist.installdir + "/mmg3d4") 

  if ( os.path.exists(mmg3d_install_dir) ):
    print "The folder %s already exists. Will not install. Exiting." % mmg3d_install_dir
    exit(1)

  os.chdir(mmg3d_src_dir)
  os.mkdir(mmg3d_build_dir)
  os.chdir(mmg3d_build_dir)

  #print "Working directory: %s" % os.getcwd()
  #print "Configuring mmg3d build ..."
  cmnd = "cmake -D CMAKE_BUILD_TYPE:STRING=Release \
 -D CMAKE_INSTALL_PREFIX:PATH=" + mmg3d_install_dir + " \ + -D COMPIL_STATIC_LIBRARY:BOOL=ON \
 -D USE_SCOTCH:BOOL=ON \
 -D INCLUDE_SCOTCH:PATH=" + optlist.installdir + "/scotch/include \
 -D LIBS_SCOTCH:PATH=" + optlist.installdir + "/scotch/lib/libscotch.a \
 -D LIBS_SCOTCHERR:PATH=" + optlist.installdir + "/scotch/lib/libscotcherr.a " + mmg3d_src_dir
  call(cmnd,shell=True)

  cmnd = "make"
  call(cmnd,shell=True)


  os.mkdir(mmg3d_install_dir)
  mmg3d_include_dir = mmg3d_install_dir + "/include"
  os.mkdir(mmg3d_include_dir)
  cmnd="cp " + mmg3d_src_dir + "/src/common/libmmg3d.h " + mmg3d_install_dir + "/include"
  call(cmnd,shell=True)
  mmg3d_lib_dir     = mmg3d_install_dir + "/lib"
  os.mkdir(mmg3d_lib_dir)
  cmnd = "cp libmmg3dlib4.0.a " + mmg3d_lib_dir
  call(cmnd,shell=True)

  os.chdir(rootpath)
  
  print " ... mmg3d installed"

# ========================================================================================================

def install_petsc(optlist):
  print "Starting installation of petsc ... "
  
  rootpath = os.getcwd()

  petsc_version = "3.2-p7"
  petsc_archive = "petsc-" + petsc_version + ".tar.gz" 

  petsc_install_dir = os.path.abspath(optlist.installdir + "/petsc-"+petsc_version)

  if ( os.path.exists(petsc_install_dir) ):
    print "The folder %s already exists. Will not install petsc. Exiting." % petsc_install_dir
    exit(1)

  os.chdir(optlist.repodir)
  cmnd="cp " + petsc_archive + " " + optlist.installdir
  call(cmnd,shell=True)
  os.chdir(optlist.installdir)
  cmnd="tar -xzf " + petsc_archive
  call(cmnd,shell=True)
  os.chdir(petsc_install_dir)

  os.environ["PETSC_ARCH"] = "linux-gnu-cxx-opt"
  os.environ["PETSC_DIR"] = petsc_install_dir

  petsc_config_opts=\
  "PETSC_DIR=" + petsc_install_dir + " PETSC_ARCH=linux-gnu-cxx-opt --COPTFLAGS=O2 --CXXOPTFLAGS=-O2 --CFLAGS=\"-O2 -fpermissive\" \
  --CXXFLAGS=\"-O2 -fpermissive\" -with-clanguage=cxx --with-c++-support --with-shared-libraries=1 --with-mpi-dir=" + optlist.depsroot + "/openmpi \
  --with-numpy --with-blas-lapack-lib-dir=" + optlist.lapack_lib_dir + " --download-parmetis=1 \
  --download-umfpack=1 --download-hypre=1 --download-mumps=1 --download-scalapack=1 \
  --download-blacs=1 --download-pastix=1 --download-ptscotch=1"

  print "PETSc config options: %s" % petsc_config_opts
 
  os.chdir(petsc_install_dir)
  cmnd="./config/configure.py " + petsc_config_opts
  call(cmnd,shell=True)


  call(optlist.makecmnd,shell=True)

  print " ... petsc installed"
  
# ========================================================================================================
  
def install_taucs(optlist):
  
# Very good information regarding compilation of Taucs can be found on: 
# http://www.resistivity.net/index.php?taucs&type=1
# http://cgal-discuss.949826.n4.nabble.com/Taucs-installation-problem-td3027923.html

# To clean build:
# rm configurator/configurator
# make clean
# rm makefile
# Configure:
#
# export OSTYPE=linux_Smoula (the name of this file without the 'mk' extension)
#
# ./configure prefix=$HOME/local/taucs
#
# make
#
# Create folders $HOME/local/taucs/include
#                $HOME/local/taucs/lib64
#
# make install

  print "Starting installation of taucs ... "
  
  rootpath = os.getcwd()
  
  taucs_archive = "taucs_full_scfd_mod.tar.gz"
  taucs_src_dir = os.path.abspath(optlist.repodir + "/taucs_full_scfd_mod")
  taucs_config_file = taucs_src_dir + "/config/linux_scfd.mk"
  
  if ( os.path.exists(optlist.taucs_install_dir) ):
    print "The folder %s already exists. Will not install taucs. Exiting." % optlist.taucs_install_dir
    exit(1)

  os.chdir(optlist.repodir)
  
  
  cmnd="tar -xzf " + taucs_archive
  call(cmnd,shell=True)
  
  os.chdir(taucs_src_dir)
  
  
    # Create new customized file UFconfig.mk
  outfile = open(taucs_config_file,'w')

  outputstr = """OBJEXT=.o
LIBEXT=.a
EXEEXT=
F2CEXT=.f
PATHSEP=/
DEFFLG=-D

FC        = gfortran
FFLAGS    = -O3 -fno-second-underscore -Wall -fPIC
FOUTFLG   = -o ./

COUTFLG   = -o ./
#CFLAGS    = -O3 -Wall -std=c99 -pedantic -lgfortran -I$(HOME)/local/cilk/include/cilk -fPIC
CFLAGS    = -O3 -Wall -std=c99 -pedantic -lgfortran -fPIC

LD        = $(CC)
LDFLAGS   =
LOUTFLG   = $(COUTFLG)

AR        = ar cr
AOUTFLG   =

RANLIB    = ranlib
RM        = rm -rf

#LIBBLAS   = -L $(HOME)/local/lapack/lib -lblas
#LIBLAPACK = -L $(HOME)/local/lapack/lib -llapack
"""

  outfile.write(outputstr)

  outputstr = "LIBBLAS   = -L " + optlist.lapack_lib_dir + " -lblas\n" + \
              "LIBLAPACK = -L " + optlist.lapack_lib_dir + " -llapack\n\n" + \
              "LIBMETIS  = -L " + optlist.parmetis_lib_dir + " -lmetis\n\n" + \
              "LIBF77 = -lgfortran -lpthread\n" + \
              "LIBC   = -lm\n"

  outfile.write(outputstr)
  
  outfile.close()

  
  
  os.chdir(taucs_src_dir)
  
  call("rm configurator/configurator",shell=True)
  call("make clean",shell=True)
  call("rm makefile",shell=True)

  os.environ["OSTYPE"] = "linux_scfd" # The name of the config file without the 'mk' extension
  call(cmnd,shell=True)
  
  cmnd = "./configure prefix=" + optlist.taucs_install_dir
  call(cmnd,shell=True)
  call("make",shell=True)
  
  os.mkdir(os.path.abspath(optlist.taucs_install_dir))
  os.mkdir(os.path.abspath(optlist.taucs_install_dir + "/include"))
  os.mkdir(os.path.abspath(optlist.taucs_install_dir + "/lib64"))
  
  call("make install",shell=True)
  
  os.chdir(rootpath)
  
  print " ... taucs installed"

# ========================================================================================================


script_usage_str = """
%s [-h | --help] [-d dir | --installdir dir] [-r | --repodir] package1 package2  ...
       installs packages package1 package2 ... placed in repository dir repoodir to installdir

----------------------------------------------------------------------------------------------------------------
Example usage: %s --repodir=$HOME/temp --installdir=$HOME/local/scfd_deps --packageset=basic
----------------------------------------------------------------------------------------------------------------
""" % (sys.argv[0],sys.argv[0])


class OptList:
  def __init__(self):
    self.repodir              = os.path.abspath("") 
    self.installdir           = os.path.abspath("")
    self.depsroot             = os.path.abspath("")
    self.makecmnd             = ""
    self.fortrancompiler      = ""
    self.ccompiler            = ""
    self.cppcompiler          = ""
    self.mpiccompiler         = ""
    self.mpicppcompiler       = ""
    self.mpi_include_dir      = os.path.abspath("")
    self.mpi_lib_dir          = os.path.abspath("")
    self.mpi_install_dir      = os.path.abspath("")
    self.cmake_install_dir    = os.path.abspath("")
    self.lapack_lib_dir       = os.path.abspath("")
    self.lapack_install_dir   = os.path.abspath("")
    self.parmetis_dir         = os.path.abspath("")
    self.parmetis_include_dir = os.path.abspath("")
    self.parmetis_lib_dir     = os.path.abspath("")
    self.parmetis_install_dir = os.path.abspath("")
    self.scotch_dir           = os.path.abspath("")
    self.scotch_include_dir   = os.path.abspath("")
    self.scotch_lib_dir       = os.path.abspath("")
    self.scotch_install_dir   = os.path.abspath("")
    self.hdf5_install_dir     = os.path.abspath("")
    self.hdf5_use_old_api     = "no"
    self.taucs_install_dir    = os.path.abspath("")

  def analyze_options(self,options):

    # REPOSITORY WITH PACKAGES:
    if options.repodir == None:
       print "Error: no repository directory was given!. Use the option --repodir"
       sys.exit(1)
    self.repodir = os.path.abspath(options.repodir) 

    # COMMON INSTALLATION DIRECTORY FOR ALL PACKAGES:
    if options.installdir == None:
       print "Error: no installation directory was given!. Use the option --installdir"
       sys.exit(1)
    self.installdir = os.path.abspath(options.installdir)
    
    if options.depsroot == None:
      print "No root repository (depsroot) given. Assuming that depsroot = installdir"
      self.depsroot = os.path.abspath(self.installdir)
    else:
      self.depsroot = os.path.abspath(self.depsroot)

    if options.makeopts == None:
      self.makecmnd = "make"
    else:
      self.makecmnd = "make " + options.makeopts

    # COMPILERS:
    if options.fortrancompiler == None:
      self.fortrancompiler = "gfortran"
    else:
      self.fortrancompiler = options.fortrancompiler

    if options.ccompiler == None:
      self.ccompiler = "gcc"
    else:
      self.ccompiler = options.ccompiler
    
    if options.cppcompiler == None:
      self.cppcompiler = "g++"
    else: 
      self.cppcompiler = options.cppcompiler
      
    if options.mpiccompiler == None:
      self.mpiccompiler = "mpicc"
    else: 
      self.mpiccompiler = options.mpiccompiler

    if options.mpicppcompiler == None:
      self.mpicppcompiler = "mpicxx"
    else: 
      self.mpicppcompiler = options.mpicppcompiler

    # PATH TO MPI HEADERS AND LIBRARIES:
    if options.mpi_include_dir == None:
      self.mpi_include_dir = os.path.abspath(options.installdir + "/openmpi/include")
    else:
      self.mpi_include_dir = os.path.abspath(options.mpi_include_dir)
 
    if options.mpi_lib_dir == None:
      self.mpi_lib_dir = os.path.abspath(options.installdir + "/openmpi/lib")
    else:
      self.mpi_lib_dir = os.path.abspath(options.mpi_lib_dir)

    # INSTALLATION DIRECTORY FOR OPENMPI:
    if options.mpi_install_dir == None:
      self.mpi_install_dir = os.path.abspath(options.installdir + "/openmpi")
    else:
      self.mpi_install_dir = os.path.abspath(options.mpi_install_dir)

    # INSTALLATION DIRECTORY FOR CMAKE:
    if options.cmake_install_dir == None:
      self.cmake_install_dir = os.path.abspath(options.installdir + "/cmake")
    else:
      self.cmake_install_dir = os.path.abspath(options.cmake_install_dir)

    # LAPACK:
    if options.lapack_lib_dir == None:
      self.lapack_lib_dir = os.path.abspath(options.installdir + "/lapack/lib")
    else:
      self.lapack_lib_dir = os.path.abspath(options.lapack_lib_dir)

    if options.lapack_install_dir == None:
      self.lapack_install_dir = os.path.abspath(options.installdir + "/lapack")
    else:
      self.lapack_install_dir = os.path.abspath(options.lapack_install_dir)
     
    # PARMETIS 
    if options.parmetis_dir == None:
      self.parmetis_dir = os.path.abspath(options.installdir + "/parmetis")
    else:
      self.parmetis_dir = os.path.abspath(options.parmetis_dir)
      
    if options.parmetis_include_dir == None:
      self.parmetis_include_dir = os.path.abspath(self.parmetis_dir + "/include")
    else:
      self.parmetis_include_dir = os.path.abspath(options.parmetis_include_dir)
      
    if options.parmetis_lib_dir == None:
      self.parmetis_lib_dir = os.path.abspath(self.parmetis_dir + "/lib")
    else:
      self.parmetis_lib_dir = os.path.abspath(options.parmetis_lib_dir)    
      
    if options.parmetis_install_dir == None:
      self.parmetis_install_dir = os.path.abspath(options.installdir + "/parmetis")
    else:
      self.parmetis_install_dir = os.path.abspath(options.parmetis_install_dir)    

    # SCOTCH:
    if options.scotch_dir == None:
      self.scotch_dir = os.path.abspath(options.installdir + "/scotch")
    else:
      self.scotch_dir = os.path.abspath(options.scotch_dir)

    if options.scotch_include_dir == None:
      self.scotch_include_dir = os.path.abspath(self.scotch_dir + "/include")
    else:
      self.scotch_include_dir = os.path.abspath(options.scotch_include_dir)
      
    if options.scotch_lib_dir == None:
      self.scotch_lib_dir = os.path.abspath(self.scotch_dir + "/lib")
    else:
      self.scotch_lib_dir = os.path.abspath(options.scotch_lib_dir)    
      
    if options.scotch_install_dir == None:
      self.scotch_install_dir = os.path.abspath(options.installdir + "/scotch")
    else:
      self.scotch_install_dir = os.path.abspath(options.scotch_install_dir)    

    # HDF:
    if options.hdf5_install_dir == None:
      self.hdf5_install_dir = os.path.abspath(options.installdir + "/hdf5")
    else:
      self.hdf5_install_dir = os.path.abspath(options.hdf5_install_dir)    

    if options.hdf5_use_old_api != None:
       self.hdf5_use_old_api = options.hdf5_use_old_api
       
    # TAUCS:
    if options.taucs_install_dir == None:
      self.taucs_install_dir = os.path.abspath(options.installdir + "/taucs")
    else:
      self.taucs_install_dir = os.path.abspath(options.taucs_install_dir)

  def print_options(self):
    print "Package repository = %s" % self.repodir
    print "Installation directory = %s" % self.installdir
    print "Root directory of all dependencies = %s" % self.depsroot
    print "Make command used for compilation of packages: %s" % self.makecmnd

# ========================================================================================================

#Langtangen, section 8.1
parser = OptionParser(usage=script_usage_str)

# help message is automatically provided, no need to specify --help or -h
parser.add_option('-r','--repodir',              dest = 'repodir',              help='repository directory with packages to compile')
parser.add_option('-i','--installdir',           dest = 'installdir',           help='installation directory')
parser.add_option('-d','--depsroot',             dest = 'depsroot',             help='root directory with other dependencies')
parser.add_option('-m','--makeopts',             dest = 'makeopts',             help='extra options for make command')
parser.add_option('-f','--fortrancompiler',      dest = 'fortrancompiler',      help='Fortran language compiler')
parser.add_option('-c','--ccompiler',            dest = 'ccompiler',            help='C language compiler')
parser.add_option('-C','--cppcompiler',          dest = 'cppcompiler',          help='C++ language compiler')
parser.add_option(     '--mpiccompiler',         dest = 'mpiccompiler',         help='MPI C language compiler')
parser.add_option(     '--mpicppcompiler',       dest = 'mpicppcompiler',       help='MPI C++ language compiler')
parser.add_option(     '--mpi-include-dir',      dest = 'mpi_include_dir',      help='Directory with MPI header files (mpi.h, ...)')
parser.add_option(     '--mpi-lib-dir',          dest = 'mpi_lib_dir',          help='Directory with MPI libraries (libmpi.a, ...)')
parser.add_option(     '--mpi-install-dir',      dest = 'mpi_install_dir',      help='Directory in which mpi should be installed')
parser.add_option(     '--cmake-install-dir',    dest = 'cmake_install_dir',    help='Installation directory for cmake')
parser.add_option(     '--lapack-lib-dir',       dest = 'lapack_lib_dir',       help='Directory containing liblapack.[so|a]')
parser.add_option(     '--lapack-install-dir',   dest = 'lapack_install_dir',   help='Installation directory for lapack')
parser.add_option(     '--parmetis-dir',         dest = 'parmetis_dir',         help='Directory containing parmetis headers and library files')
parser.add_option(     '--parmetis-include-dir', dest = 'parmetis_include_dir', help='Directory containing headers for [par]metis')
parser.add_option(     '--parmetis-lib-dir',     dest = 'parmetis_lib_dir',     help='Directory containing libraries for [par]metis')
parser.add_option(     '--parmetis-install-dir', dest = 'parmetis_install_dir', help='Installation directory for parmetis')
parser.add_option(     '--scotch-dir',           dest = 'scotch_dir',           help='Directory containing scotch headers and library files')
parser.add_option(     '--scotch-include-dir',   dest = 'scotch_include_dir',   help='Directory containing headers for [pt]scotch')
parser.add_option(     '--scotch-lib-dir',       dest = 'scotch_lib_dir',       help='Directory containing libraries for [pt]scotch')
parser.add_option(     '--scotch-install-dir',   dest = 'scotch_install_dir',   help='Installation directory for PTScotch')
parser.add_option(     '--hdf5-install-dir',     dest = 'hdf5_install_dir',     help='Installation directory for hdf5')
parser.add_option(     '--hdf5-use-old-api',     dest = 'hdf5_use_old_api',     help='Use the api of hdf5 version 1.6')
parser.add_option(     '--taucs-install-dir',    dest = 'taucs_install_dir',    help='Installation directory for Taucs')

(options, args) = parser.parse_args(sys.argv[1:])
#print options
#print args

optcontainer = OptList()
optcontainer.analyze_options(options)
#optcontainer.print_options()

if ( os.path.exists(optcontainer.installdir) ):
  print "The folder %s already exists" % optcontainer.installdir 
else:
  os.mkdir(optcontainer.installdir)

print "The repository dir is %s" % optcontainer.repodir

#install_cmake(optcontainer)
install_openmpi(optcontainer)
#install_atlas(optcontainer)
#install_suitesparse(optcontainer)
#install_lapack(optcontainer)
#install_ptscotch(optcontainer)
#install_boost(optcontainer)
#install_hdf5(optcontainer)
#install_parmetis(optcontainer)
#install_mmg3d(optcontainer)
#install_petsc(optcontainer)
#install_taucs(optcontainer)
#call(["mkdir",optcontainer.repodir])

