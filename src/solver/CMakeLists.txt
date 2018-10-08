include(PDEKITFiles.cmake)

# Collect common headers and prepend the full path
foreach(header_file ${PDEKIT_Solver_HEADERS})
 list(APPEND PDEKIT_Solver_Common_HEADERS "${PROJECT_SOURCE_DIR}/src/solver/${header_file}")
endforeach()

# Collect common headers and prepend the full path
foreach(source_file ${PDEKIT_Solver_SOURCES})
 list(APPEND PDEKIT_Solver_Common_SOURCES "${PROJECT_SOURCE_DIR}/src/solver/${source_file}")
endforeach()

add_library(solver SHARED ${PDEKIT_Solver_HEADERS} ${PDEKIT_Solver_SOURCES})
target_compile_options(solver PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(solver common mesh interpolation)

add_subdirectory( fe )
add_subdirectory( rdm )
add_subdirectory( time )
