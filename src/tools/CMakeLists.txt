#CmakeLists.txt in Tools dir

#Adds a library called Tools (libTools.a under Linux) 

include(PDEKITFiles.cmake)

ADD_LIBRARY(tools SHARED ${PDEKIT_Tools_HEADERS})

#FILE(GLOB Tools "*.hh")
# FILE(GLOB Tools "Loop.hh" "SingletonShell.hh" "FromToString.hh")

SET_TARGET_PROPERTIES(tools PROPERTIES LINKER_LANGUAGE CXX)

