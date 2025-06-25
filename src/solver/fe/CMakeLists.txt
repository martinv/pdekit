#CmakeLists.txt in RDM dir

include(PDEKITFiles.cmake)

add_library(fe SHARED ${PDEKIT_Solver_FE_HEADERS} ${PDEKIT_Solver_FE_SOURCES} )

target_link_libraries(fe mesh interpolation physics solver linear_system ${CMAKE_THREAD_LIBS_INIT})

set_target_properties(fe PROPERTIES LINKER_LANGUAGE CXX)
