# BLAS ----------------------------
if( NOT PDEKIT_HAVE_BLAS )

  if( EXISTS ${PDEKIT_BLAS_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_BLAS_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_BLAS_DIR}/lib )
    pdekit_add_trial_library_path( ${PDEKIT_BLAS_DIR}/lib64 )
  endif()

  if( EXISTS $ENV{BLAS_HOME} )
    pdekit_add_trial_library_path( $ENV{BLAS_HOME}/lib )
    pdekit_add_trial_library_path( $ENV{BLAS_HOME}/lib64 )
  endif()

  find_library(PDEKIT_BLAS_LIBRARIES blas ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(PDEKIT_BLAS_LIBRARIES blas )

  if( PDEKIT_BLAS_LIBRARIES )
    include_directories( "${PDEKIT_BLAS_DIR}/include" )
    set( PDEKIT_HAVE_BLAS 1 CACHE BOOL "Found BLAS library" )
  else()
    set( PDEKIT_HAVE_BLAS 0 )
  endif()

endif()

pdekit_log_file( "PDEKIT_HAVE_BLAS: [${PDEKIT_HAVE_BLAS}]" )
if(PDEKIT_HAVE_BLAS)
   pdekit_log_file( "  BLAS_LIBRARIES (native):     [${PDEKIT_BLAS_LIBRARIES}]" )
endif()

# LAPACK ----------------------------
if( NOT PDEKIT_HAVE_LAPACK )

  pdekit_set_trial_include_path("") # clear include search path
  pdekit_set_trial_library_path("") # clear library search path


  if( EXISTS ${PDEKIT_LAPACK_DIR} )
    pdekit_add_trial_library_path ( ${PDEKIT_LAPACK_DIR}  )
    pdekit_add_trial_library_path ( ${PDEKIT_LAPACK_DIR}/lib )
    pdekit_add_trial_library_path ( ${PDEKIT_LAPACK_DIR}/lib64 )
  endif()

  if( EXISTS $ENV{LAPACK_HOME} )
    pdekit_add_trial_library_path( $ENV{LAPACK_HOME} )
    pdekit_add_trial_library_path( $ENV{LAPACK_HOME}/lib64 )
  endif()

  
  find_library(NETLIB_LAPACK_LIBS NAMES lapack PATHS ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(NETLIB_LAPACK_LIBS NAMES lapack)

  find_library(NETLIB_LAPACKE_LIBS NAMES lapacke PATHS ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(NETLIB_LAPACKE_LIBS NAMES lapacke)


  set(PDEKIT_LAPACK_LIBRARIES "${NETLIB_LAPACK_LIBS};${NETLIB_LAPACKE_LIBS}")

endif()

if( PDEKIT_LAPACK_LIBRARIES )
  include_directories( "${PDEKIT_LAPACK_DIR}/include" )
  set( PDEKIT_HAVE_LAPACK 1 CACHE BOOL "Found LAPACK library")
else()
  set( PDEKIT_HAVE_LAPACK 0 )
endif()

if ( PDEKIT_HAVE_BLAS AND PDEKIT_HAVE_LAPACK )
  set( PDEKIT_BLAS_LAPACK_PROVIDER "Netlib" CACHE STRING "Provider of BLAS and LAPACK libraries" )
endif()

pdekit_log_file( "PDEKIT_HAVE_LAPACK: [${PDEKIT_HAVE_LAPACK}]" )
if(PDEKIT_HAVE_LAPACK)
  pdekit_log_file( "  LAPACK_LIBRARIES (native):   [${PDEKIT_LAPACK_LIBRARIES}]" )
endif()
