# ============================================================================
# This cmake file was taken from coolfluid3
# ============================================================================
# finding Boost (essential)

include( DetectBoost )

# ============================================================================
# finding MPI (essential) must be after the GlobalOptions

include( DetectMPI )
# 
# coolfluid_set_package( PACKAGE MPI DESCRIPTION "parallel communication" )

# ============================================================================
# finding Blas and Lapack (essential)

find_package(BlasLapack)

#=============================================================================

find_package(Threads) # POSIX Threads library

# ============================================================================
# find non essential packages
# ============================================================================

find_package( Trilinos )
find_package( Scotch )

# ============================================================================

# packages that don't influence functionality
# user should not be concerned with these
# so make them quiet and they don't appear on the enable package list

find_package(OpenSSL           QUIET ) # OpenSSL library
find_package(PythonInterp      QUIET ) # Python interpreter
find_package(CMath             QUIET ) # find the math library
find_package(Realtime          QUIET ) # POSIX Realtime library
find_package(Valgrind          QUIET ) # valgrind for profiling and memmory leak detection
find_package(GooglePerftools   QUIET ) # dynamic profiler and memory checker
#find_package(Threads           QUIET ) # POSIX Threads library
#find_package(PThread           QUIET ) # POSIX Threads library
