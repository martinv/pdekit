#CmakeLists.txt in RDM dir

#Adds a library called RDM (libRDM.a under Linux) from the source files

include(PDEKITFiles.cmake)

# Collect common headers and prepend the full path
#foreach(header_file ${PDEKIT_Solver_HEADERS})
# list(APPEND PDEKIT_Solver_Common_HEADERS "${PROJECT_SOURCE_DIR}/src/solver/${header_file}")
#endforeach()

#add_library(rdm SHARED ${PDEKIT_Solver_Common_HEADERS} ${PDEKIT_Solver_Common_SOURCES}
#                       ${PDEKIT_Solver_RDM_HEADERS}    ${PDEKIT_Solver_RDM_SOURCES})

add_library(rdm SHARED ${PDEKIT_Solver_RDM_HEADERS} ${PDEKIT_Solver_RDM_SOURCES})

target_compile_options(rdm PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
# Note that 'solver' is an independent library that is linked to RDM
target_link_libraries(rdm mesh interpolation physics solver linear_system ${CMAKE_THREAD_LIBS_INIT})

#set_target_properties(rdm PROPERTIES LINKER_LANGUAGE CXX)

add_custom_command(TARGET rdm
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/src/solver/rdm/RDSolverOptions.dat  ${PROJECT_BINARY_DIR}/examples/testcases-rdm)
#                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/src/solver/rdm/RDSolverOptions.dat  ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR})
