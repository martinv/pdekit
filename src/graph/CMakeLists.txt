#CmakeLists.txt in graph dir

include(PDEKITFiles.cmake)

add_library(graph SHARED ${PDEKIT_Graph_HEADERS} ${PDEKIT_Graph_SOURCES})

target_compile_options(graph PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(graph ${PDEKIT_SCOTCH_LIBRARIES} pthread)

install(TARGETS graph DESTINATION lib)

