#CmakeLists.txt in common dir

# Adds a library called common (libcommon.a under Linux)
# from the source files listed below

include(PDEKITFiles.cmake)
include(MPI/PDEKITFiles.cmake)

add_library(common SHARED ${PDEKIT_Common_HEADERS}
                          ${PDEKIT_Common_SOURCES}
                          ${PDEKIT_MPI_HEADERS}
                          ${PDEKIT_MPI_SOURCES})

target_compile_options(common PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(common ${MPI_CXX_LIBRARIES})

install(FILES ${PDEKIT_Common_HEADERS} DESTINATION include/common)
install(FILES ${PDEKIT_MPI_HEADERS} DESTINATION include/common)
install(TARGETS common DESTINATION lib)


