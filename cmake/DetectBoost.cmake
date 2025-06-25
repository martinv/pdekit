# =============================================================================================================
# This cmake file was taken from coolfluid3
# Look for boost. If PDEKIT_DEPS_ROOT is defined, try to link to boost located inside PDEKIT_DEPS_ROOT
# This is done via appending PDEKIT_DEPS_ROOT to CMAKE_PREFIX_PATH in ResetVariables.cmake
# =============================================================================================================

set( Boost_USE_STATIC_LIBS ${CF_ENABLE_STATIC} )
set( Boost_USE_MULTITHREAD ON  )
# find based on minimal version defined below
set( Boost_FIND_VERSION        ON   )
set( Boost_FIND_VERSION_MAJOR  "1"  )
set( Boost_FIND_VERSION_MINOR  "59" )
set( Boost_FIND_VERSION_PATCH  "0"  )
# older cmakes dont have these versions
set( Boost_ADDITIONAL_VERSIONS "1.59" "1.59.0" "1.58" "1.58.0" "1.57" "1.57.0" "1.56" "1.56.0" "1.55" "1.55.0" "1.54" "1.54.0" "1.53" "1.53.0" "1.52" "1.52.0" )
# components to search for

list( APPEND PDEKIT_Boost_COMPONENTS thread iostreams filesystem system regex unit_test_framework date_time program_options )

find_package( Boost COMPONENTS ${PDEKIT_Boost_COMPONENTS} QUIET )

#add_definitions( -DBOOST_ENABLE_ASSERT_HANDLER )

# if not found give more information
if( NOT Boost_FOUND )
  pdekit_log( ${Boost_ERROR_REASON} )
  message( FATAL_ERROR "Boost is required to compile PDEKIT" )
endif()

# add boost include path and link path
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories(${Boost_LIBRARY_DIRS})

# add boost libraries to list of third party libraries
list( APPEND PDEKIT_DEPS_LIBRARIES ${Boost_LIBRARIES} )

# filter out the unit test libs from the boost libraries
# only unit tests link to this
set(PDEKIT_BOOST_LIBRARIES "" )
foreach( blib ${Boost_LIBRARIES} )
  if( NOT ${blib} MATCHES "[a-zA-Z0-9]*unit_test_framework[a-zA-Z0-9]*" )
    list( APPEND PDEKIT_BOOST_LIBRARIES ${blib} )
  endif()
endforeach()


message("=== Boost include dirs: ${Boost_INCLUDE_DIRS}")
message("=== Boost library dirs: ${Boost_LIBRARY_DIRS}")
