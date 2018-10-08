if( NOT PDEKIT_HAVE_BLAS )

  if( EXISTS ${PDEKIT_BLAS_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_BLAS_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_BLAS_DIR}/lib )
  endif()

  if( EXISTS $ENV{BLAS_HOME} )
    pdekit_add_trial_library_path( $ENV{BLAS_HOME}/lib )
  endif()

  find_library(PDEKIT_BLAS_LIBRARIES openblas ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(PDEKIT_BLAS_LIBRARIES openblas )

  unset(PDEKIT_LAPACK_LIBRARIES)

  if( PDEKIT_BLAS_LIBRARIES )
    include_directories( ${PDEKIT_BLAS_DIR}/include )
    set( PDEKIT_HAVE_BLAS 1   CACHE BOOL "Found BLAS library" )
    set( PDEKIT_HAVE_LAPACK 1 CACHE BOOL "Found LAPACK library" )
  else()
    set( PDEKIT_HAVE_BLAS 0 )
    set( PDEKIT_HAVE_LAPACK 0 )
  endif()

endif()

pdekit_log_file( "PDEKIT_HAVE_BLAS:   [${PDEKIT_HAVE_BLAS}]" )
if(PDEKIT_HAVE_BLAS)
   pdekit_log_file( "  BLAS_LIBRARIES:     [${PDEKIT_BLAS_LIBRARIES}]" )
endif()
pdekit_log_file( "PDEKIT_HAVE_LAPACK: [${PDEKIT_HAVE_LAPACK}]" )
if(PDEKIT_HAVE_LAPACK)
   pdekit_log_file( "  LAPACK_LIBRARIES:     [${PDEKIT_LAPACK_LIBRARIES}]" )
endif()
