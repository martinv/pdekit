# The following variables will be defined at the end of this output:
# PDEKIT_TRILINOS_LIBS
# PDEKIT_TRILINOS_INCLUDE_DIRS
# PDEKIT_HAVE_TRILINOS

if ( PDEKIT_USE_TRILINOS )
  if ( DEFINED PDEKIT_TRILINOS_DIR AND NOT ${PDEKIT_TRILINOS_DIR} STREQUAL "") # Here we should not dereference the variable SCFD_TRILINOS_DIR !!!
    message("Trilinos prefix = ${PDEKIT_TRILINOS_DIR}")

    find_package( Trilinos PATHS ${PDEKIT_TRILINOS_DIR}/lib/cmake/Trilinos ${PDEKIT_TRILINOS_DIR}/include)

    if( Trilinos_FOUND )

      message("\nFound Trilinos!  Here are the details: ")
      message("   Trilinos_DIR = ${Trilinos_DIR}")
      message("   Trilinos_VERSION = ${Trilinos_VERSION}")
      #message("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
      #message("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
      #message("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
      #message("   Trilinos_LIBRARY_DIRS = ${Trilinos_LIBRARY_DIRS}")
      #message("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
      #message("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
      #message("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
      #message("   Trilinos_TPL_LIBRARY_DIRS = ${Trilinos_TPL_LIBRARY_DIRS}")
      message("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
      message("End of Trilinos details\n")

      set(PDEKIT_TRILINOS_LIBS "")

      foreach(lib_name ${Trilinos_LIBRARIES})
        #message("TESTLIB = ${lib_name}")
        # The found library will be stored in  the variable ${lib_name}_lib
        # The name of the variable is different for each detected library.
        # This is necessary, as find_library caches the variable in which it stores
        # found libraries and using the same variable in loop
        # results in finding the same library over and over

        find_library( ${lib_name}_lib ${lib_name} PATHS  ${Trilinos_LIBRARY_DIRS}  NO_DEFAULT_PATH )

        #message(" --> ${${lib_name}_lib}")
        list( APPEND PDEKIT_TRILINOS_LIBS ${${lib_name}_lib} )
        mark_as_advanced( ${lib_name}_lib )

      endforeach()

      list(APPEND PDEKIT_TRILINOS_LIBS ${Trilinos_TPL_LIBRARIES} )

      #message("Printing the list of trilinos libs:")
      #message( "${PDEKIT_TRILINOS_LIBS}" )

      set( PDEKIT_TRILINOS_INCLUDE_DIRS "" )
      list( APPEND PDEKIT_TRILINOS_INCLUDE_DIRS ${Trilinos_INCLUDE_DIRS} )
      list( APPEND PDEKIT_TRILINOS_INCLUDE_DIRS ${Trilinos_TPL_INCLUDE_DIRS} )
      #message("PDEKIT Trilinos include dirs = ${PDEKIT_TRILINOS_INCLUDE_DIRS}")


      include_directories( ${PDEKIT_TRILINOS_INCLUDE_DIRS} )

      set( PDEKIT_HAVE_TRILINOS 1 )

    endif( Trilinos_FOUND )


  else()

    message(FATAL_ERROR "To search for trilinos, the variable PDEKIT_TRILINOS_DIR has to be set!")

  endif( DEFINED PDEKIT_TRILINOS_DIR AND NOT ${PDEKIT_TRILINOS_DIR} STREQUAL "" )

else()

  set( PDEKIT_HAVE_TRILINOS 0 )

endif( PDEKIT_USE_TRILINOS )
