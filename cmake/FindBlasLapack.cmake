# =========================================
# This cmake file was taken from coolfluid3
# =========================================

# Confirm that liblapack library is installed
# This module defines
#   PDEKIT_HAVE_LAPACK
#   PDEKIT_HAVE_BLAS
#   PDEKIT_HAVE_BLAS_LAPACK
#   BLAS_LAPACK_LIBRARIES
# Require both lapack and blas

# Get the list of recognized Blas and Lapack providers
# This property is set in DefineGlobalOptions.cmake
get_property(AVAILABLE_BLAS_LAPACK_PROVIDERS CACHE PDEKIT_BLAS_LAPACK_PROVIDER PROPERTY STRINGS)

# Check if the configured Blas/Lapack provider is a valid one
check_if_arg_value_is_valid("${AVAILABLE_BLAS_LAPACK_PROVIDERS}" ${PDEKIT_BLAS_LAPACK_PROVIDER} result)

if ( NOT ${result} )
  message(FATAL_ERROR "PDEKIT_BLAS_LAPACK_PROVIDER cannot be set to ${PDEKIT_BLAS_LAPACK_PROVIDER}. Use one of the following: ${AVAILABLE_BLAS_LAPACK_PROVIDERS}")
else()
  message(STATUS "Available providers for blas and lapack: ${AVAILABLE_BLAS_LAPACK_PROVIDERS}")
  message(STATUS "Configured provider: ${PDEKIT_BLAS_LAPACK_PROVIDER}")
endif()

if( NOT PDEKIT_LAPACK_LIBRARIES ) # This variable is set by cmake using a standard search for
                                # Blas and Lapack
  # ===========================================================================
  # 1) Compile the OpenBLAS library, which should provide Blas, CBlas and
  #    Lapack functionality
  # ===========================================================================

  if ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Internal" )
    include ( BuildInternalBlasLapack )

  # ===========================================================================
  # 2) Use the default implemenation from Netlib
  # ===========================================================================

  elseif ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Netlib" )

    include ( FindNativeBlasLapack )

  # ===========================================================================
  # 3) Use OpenBLAS
  #    In this case, one library (libopenblas) contains both
  #    Blas and Lapack functionality ...
  # ===========================================================================

  elseif ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "OpenBLAS" )

    include ( FindOpenBlas )

  # ===========================================================================
  # 4) Use IntelMKL
  # ===========================================================================

  elseif ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "IntelMKL" )

    include ( FindIntelMKL )

  endif ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Internal" )


# BOTH ###########################

  if( PDEKIT_HAVE_BLAS AND PDEKIT_HAVE_LAPACK )
    set( BLAS_LAPACK_LIBRARIES   "${PDEKIT_LAPACK_LIBRARIES};${PDEKIT_BLAS_LIBRARIES}" CACHE STRING "BLAS and LAPACK libraries")
    #message("Debug: setting blas and lapack libraries to: ${BLAS_LAPACK_LIBRARIES}")
    set( PDEKIT_HAVE_BLAS_LAPACK 1 CACHE BOOL "Found BLAS and LAPACK libraries")

    #pdekit_log_file( "--> PDEKIT_HAVE_BLAS_LAPACK: [${PDEKIT_HAVE_BLAS_LAPACK}]" )
    #pdekit_log_file( "--> BLAS_LAPACK_LIBRARIES: [${BLAS_LAPACK_LIBRARIES}]" )

  endif()

  mark_as_advanced( PDEKIT_LAPACK_LIBRARIES PDEKIT_BLAS_LIBRARIES )

#################################

else()

  # user provided directly the libraries of LAPACK
  # TODO: test  that they actually work

  set(PDEKIT_HAVE_LAPACK       ON CACHE BOOL "Found LAPACK library")
  set(PDEKIT_HAVE_BLAS         ON CACHE BOOL "Found BLAS   library")
  set(PDEKIT_HAVE_BLAS_LAPACK  ON CACHE BOOL "Found BLAS and LAPACK libraries")

  set( BLAS_LAPACK_LIBRARIES   "${PDEKIT_LAPACK_LIBRARIES}" CACHE STRING "BLAS and LAPACK libraries")

  mark_as_advanced( BLAS_LAPACK_LIBRARIES PDEKIT_LAPACK_LIBRARIES )

  #pdekit_log_file( "-->> PDEKIT_HAVE_BLAS_LAPACK: [${PDEKIT_HAVE_BLAS_LAPACK}]" )
  #pdekit_log_file( "-->> BLAS_LAPACK_LIBRARIES: [${BLAS_LAPACK_LIBRARIES}]" )

endif()

mark_as_advanced( PDEKIT_HAVE_LAPACK PDEKIT_HAVE_BLAS PDEKIT_HAVE_BLAS_LAPACK BLAS_LAPACK_LIBRARIES )

#if ( ${PDEKIT_HAVE_BLAS_LAPACK} )
#    list( APPEND PDEKIT_DEPS_LIBRARIES ${BLAS_LAPACK_LIBRARIES} )
#endif()

set( BlasLapack_FOUND ${PDEKIT_HAVE_BLAS_LAPACK} )
#pdekit_set_package( PACKAGE BlasLapack 
#                       DESCRIPTION "linear algebra"
#                       VARS BLAS_LAPACK_LIBRARIES)

# MAKE SURE THE FOLLOWING VARIABLES ARE DEFINED: IF THEY DON'T EXIST,
# SET THEM TO ZERO

set( PDEKIT_BLAS_LAPACK_PROVIDER_IS_NETLIB 0 CACHE BOOL "The Blas and Lapack libs are provided by Netlib" )
set ( PDEKIT_BLAS_LAPACK_PROVIDER_IS_OPENBLAS 0 CACHE BOOL "The Blas and Lapack libs are provided by OpenBlas" )

if ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Netlib" )
  set( PDEKIT_BLAS_LAPACK_PROVIDER_IS_NETLIB 1 CACHE BOOL "The Blas and Lapack libs are provided by Netlib" FORCE )
endif( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Netlib" )

if ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "OpenBLAS" )
  set ( PDEKIT_BLAS_LAPACK_PROVIDER_IS_OPENBLAS 1 CACHE BOOL "The Blas and Lapack libs are provided by OpenBlas" FORCE )
endif( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "OpenBLAS" )

#if ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Netlib" )
#  set( PDEKIT_BLAS_LAPACK_PROVIDER_IS_NETLIB 1 CACHE BOOL "The Blas and Lapack libs are provided by Netlib" )
#else ()
#  set( PDEKIT_BLAS_LAPACK_PROVIDER_IS_NETLIB 0 CACHE BOOL "The Blas and Lapack libs are provided by Netlib" )
#endif( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "Netlib" )

#if ( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "OpenBLAS" )
#  set ( PDEKIT_BLAS_LAPACK_PROVIDER_IS_OPENBLAS 1 CACHE BOOL "The Blas and Lapack libs are provided by OpenBlas" )
#else ()
#  set( PDEKIT_BLAS_LAPACK_PROVIDER_IS_OPENBLAS 0 CACHE BOOL "The Blas and Lapack libs are provided by OpenBlas" )
#endif( ${PDEKIT_BLAS_LAPACK_PROVIDER} STREQUAL "OpenBLAS" )

