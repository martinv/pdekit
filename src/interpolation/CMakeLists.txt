#CmakeLists.txt in interpolation dir

#Adds a library called interpolation (libinterpolation.{a.so} under Linux) from the source files

include(PDEKITFiles.cmake)

add_library(interpolation SHARED ${PDEKIT_Interpolation_HEADERS} 
                                 ${PDEKIT_Interpolation_SOURCES} )

target_compile_options(interpolation PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(interpolation mesh common linear_system)

install(FILES ${PDEKIT_Interpolation_HEADERS} DESTINATION include/interpolation)
install(TARGETS interpolation DESTINATION lib)
