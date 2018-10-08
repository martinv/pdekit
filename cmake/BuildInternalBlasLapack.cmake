message(STATUS "Including internal BLAS and LAPACK")

include(ExternalProject)

ExternalProject_Add(OpenBLASInternal
  PREFIX ${CMAKE_BINARY_DIR}/contrib/OpenBLAS
  URL ${CMAKE_SOURCE_DIR}/contrib/OpenBLAS.git
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
  BUILD_IN_SOURCE 1
  INSTALL_COMMAND make install PREFIX=${CMAKE_BINARY_DIR}/contrib/OpenBLAS
)

#ExternalProject_Get_Property(OpenBLAS binary_dir)
#set(OpenBLAS_internal_build_dir ${binary_dir})
# message("The build dir for OpenBlas is ${binary_dir}")

set(OpenBLAS_internal_lib ${CMAKE_BINARY_DIR}/contrib/OpenBLAS/lib/libopenblas.so)
#find_library( OpenBLAS_internal_lib NAMES openblas PATHS ${CMAKE_BINARY_DIR}/contrib/OpenBLAS/lib NO_DEFAULT_PATH )

include_directories( ${CMAKE_BINARY_DIR}/contrib/OpenBLAS/include )

set(PDEKIT_HAVE_BLAS         ON CACHE BOOL "Found BLAS   library")
set(PDEKIT_HAVE_LAPACK       ON CACHE BOOL "Found LAPACK library")
set(PDEKIT_HAVE_BLASLAPACK   ON CACHE BOOL "Found BLAS and LAPACK libraries")

set( BLASLAPACK_LIBRARIES ${OpenBLAS_internal_lib} CACHE STRING "BLAS and LAPACK libraries")
set( PDEKIT_BLASLAPACK_PROVIDER "OpenBLASInternal" CACHE STRING "Provider of BLAS and LAPACK libraries")
