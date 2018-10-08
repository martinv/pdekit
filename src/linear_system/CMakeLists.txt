#CmakeLists.txt in linear_system dir

include(PDEKITFiles.cmake)

add_library(linear_system SHARED ${PDEKIT_Linear_System_HEADERS}
                                 ${PDEKIT_Linear_System_SOURCES})

target_compile_options(linear_system PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(linear_system common mesh ${PDEKIT_TRILINOS_LIBS} )

install(FILES ${PDEKIT_Linear_System_HEADERS} DESTINATION include/linear_system)
install(TARGETS linear_system DESTINATION lib)

