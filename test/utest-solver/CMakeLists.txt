# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

if ( PDEKIT_HAVE_TRILINOS )

set( no_mpi_proc 1 )

#set(utest-trilinos-crs-matrix_parameters -np ${no_mpi_proc} "./utest-trilinos-crs-matrix")
#add_executable(utest-trilinos-crs-matrix utest-trilinos-crs-matrix.cpp)
#add_test(NAME utest-trilinos-crs-matrix COMMAND "mpirun" ${utest-trilinos-crs-matrix_parameters} )
#target_link_libraries(utest-trilinos-crs-matrix common linear_system boost_unit_test_framework)

set(utest-L2-global-projection_parameters -np ${no_mpi_proc} "./utest-L2-global-projection")
add_executable(utest-L2-global-projection utest-L2-global-projection.cpp)
add_test(NAME utest-L2-global-projection COMMAND "mpirun" ${utest-L2-global-projection_parameters} )
target_link_libraries(utest-L2-global-projection common mesh linear_system solver boost_unit_test_framework)

endif ( PDEKIT_HAVE_TRILINOS )

add_executable(utest-coarse-scale-correction utest-coarse-scale-correction.cpp)
add_test(utest-coarse-scale-correction utest-coarse-scale-correction)
target_link_libraries(utest-coarse-scale-correction common mesh interpolation solver linear_system boost_unit_test_framework)

add_executable(utest-L2-local-projection utest-L2-local-projection.cpp)
add_test(utest-L2-local-projection utest-L2-local-projection)
target_link_libraries(utest-L2-local-projection common mesh solver boost_unit_test_framework)

add_executable(utest-residual-restriction-local utest-residual-restriction-local.cpp)
add_test(utest-residual-restriction-local utest-residual-restriction-local)
target_link_libraries(utest-residual-restriction-local common mesh solver boost_unit_test_framework)

add_executable(utest-rdm-time-update utest-rdm-time-update.cpp)
add_test(utest-rdm-time-update utest-rdm-time-update)
target_link_libraries(utest-rdm-time-update common mesh solver rdm boost_unit_test_framework)

add_executable(utest-solver-setup-algorithm  utest-solver-setup-algorithm.cpp)
add_test(utest-solver-setup-algorithm utest-solver-setup-algorithm)
target_link_libraries(utest-solver-setup-algorithm common mesh solver boost_unit_test_framework)

add_executable(utest-elemental-matrix-operator utest-elemental-matrix-operator.cpp)
add_test(utest-elemental-matrix-operator utest-elemental-matrix-operator)
target_link_libraries(utest-elemental-matrix-operator common mesh solver boost_unit_test_framework)

add_executable(utest-Euler-continuous-rdm utest-Euler-continuous-rdm.cpp)
add_test(utest-Euler-continuous-rdm utest-Euler-continuous-rdm)
target_link_libraries(utest-Euler-continuous-rdm common mesh solver rdm boost_unit_test_framework)

add_executable(utest-RD-solver utest-RD-solver.cpp)
add_test(utest-RD-solver utest-RD-solver)
target_link_libraries(utest-RD-solver common mesh solver rdm boost_unit_test_framework)

add_executable(utest-quad-adapt-consistency utest-quad-adapt-consistency.cpp)
add_test(utest-quad-adapt-consistency utest-quad-adapt-consistency)
target_link_libraries(utest-quad-adapt-consistency common mesh interpolation tools boost_unit_test_framework)

add_executable(utest-mesh-function-snapshot utest-mesh-function-snapshot.cpp)
add_test(utest-mesh-function-snapshot utest-mesh-function-snapshot)
target_link_libraries(utest-mesh-function-snapshot common mesh interpolation tools boost_unit_test_framework)

# -----------------------------------------------------------------------------

add_subdirectory(fe)

# -----------------------------------------------------------------------------

add_custom_command(TARGET utest-coarse-scale-correction
                  POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p2.msh
                                                                ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR} )

add_custom_command(TARGET utest-quad-adapt-consistency
                  POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mini.msh
                                                                ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mini_mini.msh
                                                                ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p1_tri_no_bdry.msh
                                                                ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR} )

