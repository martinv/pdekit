#=============================================================================
# This cmake file was taken from coolfluid3
#=============================================================================

##############################################################################
# include cmake macros
##############################################################################

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(CheckSymbolExists)
include(CheckFunctionExists)
include(CheckLibraryExists)
include(CheckTypeSize)
include(CheckCSourceCompiles)
include(CheckCXXSourceCompiles)
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

include(FeatureSummary)

##############################################################################
# include pdekit functions
##############################################################################

include(functions/CheckIfArgValueIsValid)


##############################################################################
# include pdekit macros
##############################################################################

include(macros/PDEKitLogToFile)
include(macros/PDEKitSearchPaths)
