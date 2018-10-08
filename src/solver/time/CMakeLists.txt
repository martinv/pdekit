#CmakeLists.txt in time dir

include(PDEKITFiles.cmake)

add_library(time SHARED ${PDEKIT_Solver_TIME_HEADERS} ${PDEKIT_Solver_TIME_SOURCES})

target_link_libraries(time linear_system)

set_target_properties(time PROPERTIES LINKER_LANGUAGE CXX)
