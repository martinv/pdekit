cmake_minimum_required( VERSION 3.5 FATAL_ERROR )

project( PDEKIT C CXX Fortran )

# =============================================================================
# Disallow in-source build
# =============================================================================
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "${PROJECT_NAME} requires an out of source build.\nPlease create a separate build directory and run 'cmake path/to/project [options]' there.")
endif()

set( CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} )

# =============================================================================
message( STATUS "[configuring ${PROJECT_NAME}]")
# =============================================================================

include( DefineMacros ) # define your own macros
include( DefineInstallPaths ) # define where to install the files after building, set the log file

# =============================================================================
# Detect operating system
# =============================================================================

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    # Linux specific code
    set(OperatingSystem "Linux")
    message ("The operating system is Linux")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# =============================================================================
# Set CMake standard and compiler flags
# =============================================================================

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  message("Using Clang")
  set(PDEKIT_CMAKE_CXX_FLAGS "-fPIC;-Wall;-pedantic-errors;-Wno-return-type-c-linkage;-Wno-unused-local-typedef;-Wno-undefined-var-template")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message("Using GCC")
  set(PDEKIT_CMAKE_CXX_FLAGS "-fPIC;-Wall;-pedantic-errors;-Wno-unused-local-typedefs;-Wno-deprecated-declarations")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message("using Intel C++")
  set(PDEKIT_CMAKE_CXX_FLAGS "-fPIC;-Wall;-pedantic-errors")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  message("Using Visual Studio C++")
endif()

# =============================================================================
# Reset some variables
# =============================================================================
include( ResetVariables )

# =============================================================================
pdekit_log(" setting build options ")
# =============================================================================

include( DefineGlobalOptions )   # user-defined options
include( CheckOperatingSystem )  # check for operating system features
include( CheckCompilerFeatures ) # check compiler features

# =============================================================================
pdekit_log(" generic configuration ")
# =============================================================================

include( DetectSoftware )

set(CMAKE_BUILD_TYPE RelWithDebInfo)


set( PDEKIT_MESH_REPO_DIR ${PDEKIT_SOURCE_DIR}/mesh_repository )

include_directories( ${PDEKIT_LAPACK_DIR}/include )
include_directories( ${PDEKIT_SOURCE_DIR} )
include_directories( ${PDEKIT_SOURCE_DIR}/src )

#enable_testing()

# configure a header file to pass some of the CMake settings
# to the source code

configure_file (
  "${PROJECT_SOURCE_DIR}/PDEKit_Config.hpp.in"
  "${PROJECT_BINARY_DIR}/PDEKit_Config.hpp"
  )

# add the binary tree to the search path for include files
# so that we will find PDEKit_Config.hpp
include_directories("${PROJECT_BINARY_DIR}")

# =============================================================================
# Subdirectories
# =============================================================================

add_subdirectory( src )
add_subdirectory( test )
add_subdirectory( examples )
