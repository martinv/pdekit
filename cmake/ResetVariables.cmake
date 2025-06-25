# reset the list of kernel libs

#set( CF_KERNEL_LIBS "" CACHE INTERNAL "" )
#set( CF_PLUGIN_LIST "" CACHE INTERNAL "" )

# reset the list of project n orphan files

#set( CF_PROJECT_FILES "" CACHE INTERNAL "" )
#set( CF_ORPHAN_FILES  "" CACHE INTERNAL "" )

# user define that affects many search paths simultaneously

if( DEFINED PDEKIT_DEPS_ROOT )
  list(APPEND CMAKE_PREFIX_PATH ${PDEKIT_DEPS_ROOT})
endif()

# finding boost with BOOST_ROOT avoids conflicts with the system libraries

if( NOT DEFINED BOOST_ROOT )
  set( BOOST_ROOT  "${PDEKIT_DEPS_ROOT}/boost" CACHE INTERNAL "" )
endif()

# finding MPI with MPI_HOME or MPI_ROOT avoids conflicts with the system libraries

if( NOT DEFINED MPI_ROOT )
  set( MPI_ROOT  ${PDEKIT_DEPS_ROOT}/openmpi CACHE INTERNAL "" )
  set( MPI_HOME  ${PDEKIT_DEPS_ROOT}/openmpi CACHE INTERNAL "" )
endif()

# temporary

#if( NOT DEFINED CF_TMP_HAVE_SIMPLECOMM )
#  option( CF_TMP_HAVE_SIMPLECOMM "Use SimpleComm" OFF )
#endif()
