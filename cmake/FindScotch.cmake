# Pastix requires SCOTCH or METIS (partitioning and reordering tools)

if ( NOT PDEKIT_SCOTCH_LIBRARIES )
  pdekit_set_trial_include_path("") # clear include search path
  pdekit_set_trial_library_path("") # clear library search path

  if ( EXISTS ${PDEKIT_SCOTCH_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_SCOTCH_DIR} )
    pdekit_add_trial_library_path( ${PDEKIT_SCOTCH_DIR}/lib )
    message("Have trial library paths: ${TRIAL_LIBRARY_PATHS}")
  endif()

  if ( EXISTS $ENV{SCOTCH_HOME} )
    pdekit_add_trial_library_path( $ENV{SCOTCH_HOME}/lib )
  endif()

  find_library(SCOTCH_LIBRARY       scotch      ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(PTSCOTCH_LIBRARY     ptscotch    ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(SCOTCHERR_LIBRARY    scotcherr   ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)
  find_library(PTSCOTCHERR_LIBRARY  ptscotcherr ${TRIAL_LIBRARY_PATHS} NO_DEFAULT_PATH)

  set(DETECTED_PDEKIT_SCOTCH_LIBRARIES
       "${SCOTCH_LIBRARY};${PTSCOTCH_LIBRARY};${SCOTCHERR_LIBRARY};${PTSCOTCHERR_LIBRARY}"
  )
#endif()

  message("Detected scotch libraries: ${DETECTED_PDEKIT_SCOTCH_LIBRARIES}")

  if ( DETECTED_PDEKIT_SCOTCH_LIBRARIES )
    include_directories( "${PDEKIT_SCOTCH_DIR}/include" )
    set( PDEKIT_HAVE_SCOTCH 1 CACHE STRING "Found SCOTCH library" )
    set( PDEKIT_SCOTCH_LIBRARIES "${DETECTED_PDEKIT_SCOTCH_LIBRARIES}" CACHE STRING "Scotch libraries" )
  else()
    #set( PDEKIT_HAVE_SCOTCH 0 CACHE BOOL "Found SCOTCH library" )
    set( PDEKIT_HAVE_SCOTCH 0 )
  endif()

endif(NOT PDEKIT_SCOTCH_LIBRARIES)

pdekit_log_file( "PDEKIT_HAVE_SCOTCH: [${PDEKIT_HAVE_SCOTCH}]" )
if(PDEKIT_HAVE_SCOTCH)
   pdekit_log_file( "  SCOTCH_LIBRARIES:     [${PDEKIT_SCOTCH_LIBRARIES}]" )
   include_directories( "${PDEKIT_SCOTCH_DIR}/include" )
endif()


#include(FindPackageHandleStandardArgs)
#find_package_handle_standard_args(
#    SCOTCH DEFAULT_MSG SCOTCH_INCLUDES SCOTCH_LIBRARIES)

#mark_as_advanced(SCOTCH_INCLUDES SCOTCH_LIBRARIES)
