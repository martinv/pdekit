###############################################################################
### Check for Linux
if( UNIX AND NOT APPLE )
  if( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
    set( PDEKIT_OS_LINUX 1 )
  else()
    set( PDEKIT_OS_UNRECOGNIZED_REASON "Unrecognized UNIX type : pdekit has only been tested for Linux or MacOSX type UNIX'es")
  endif()
endif()

###############################################################################
### Check for Apple MacOSX
#if( APPLE )
#  if( UNIX AND ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
#    set( PDEKIT_OS_MACOSX 1 )
#  else()
#    set( PDEKIT_OS_UNRECOGNIZED_REASON "Unrecognized APPLE type : pdekit has only been tested  only with Apple MacOSX ( Darwin ) systems.")
#  endif()
#endif()

###############################################################################
### Check for Windows
#if( WIN32 )
#  if( MSVC OR MINGW )
#    set( PDEKIT_OS_WINDOWS 1 )
#  else()
#    set( PDEKIT_OS_UNRECOGNIZED_REASON "Unrecognized WINDOWS type : pdekit has only been tested with Win32 and MSVC or MingGW compiler.")
#  endif()
#endif()

###############################################################################
### FINAL MESSAGE
if( PDEKIT_OS_UNRECOGNIZED_REASON )
  set( PDEKIT_OS_UNRECOGNIZED 1 )
  if( NOT PDEKIT_SKIP_OS_TEST )
    set( FULL_MSG "${PDEKIT_OS_UNRECOGNIZED_REASON} Set CMake variable PDEKIT_SKIP_OS_TEST to avoid this error" )
    message( FATAL_ERROR ${FULL_MSG} )
  else()
    message( STATUS ${PDEKIT_OS_UNRECOGNIZED_REASON} )
    message( STATUS "Nevertheless we try to continue ..." )
  endif()
endif()

###############################################################################

