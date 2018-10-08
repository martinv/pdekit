# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

set( no_mpi_proc 1 )

set(utest-trilinos-crs-matrix_parameters -np ${no_mpi_proc} "./utest-trilinos-crs-matrix")
add_executable(utest-trilinos-crs-matrix utest-trilinos-crs-matrix.cpp)
add_test(NAME utest-trilinos-crs-matrix COMMAND "mpirun" ${utest-trilinos-crs-matrix_parameters} )
target_link_libraries(utest-trilinos-crs-matrix common linear_system boost_unit_test_framework)

set(utest-trilinos-crs-mat-block-view-params -np ${no_mpi_proc} "./utest-trilinos-crs-mat-block-view")
add_executable(utest-trilinos-crs-mat-block-view utest-trilinos-crs-mat-block-view.cpp)
add_test(NAME utest-trilinos-crs-mat-block-view COMMAND "mpirun" ${utest-trilinos-crs-mat-block-view-params} )
target_link_libraries(utest-trilinos-crs-mat-block-view common linear_system boost_unit_test_framework)

set(utest-ls-trilinos_parameters -np ${no_mpi_proc} "./utest-ls-trilinos")
add_executable(utest-ls-trilinos utest-ls-trilinos.cpp)
add_test(NAME utest-ls-trilinos COMMAND "mpirun" ${utest-ls-trilinos_parameters} )
target_link_libraries(utest-ls-trilinos common linear_system boost_unit_test_framework)

if ( PDEKIT_HAVE_TRILINOS )

set(utest-tpetra-block-preconditioner-params -np ${no_mpi_proc} "./utest-tpetra-block-preconditioner")
add_executable(utest-tpetra-block-preconditioner utest-tpetra-block-preconditioner)
add_test(NAME utest-tpetra-block-preconditioner COMMAND "mpirun" ${utest-tpetra-block-preconditioner-params} )
target_link_libraries(utest-tpetra-block-preconditioner common linear_system boost_unit_test_framework)

set(utest-Poisson_parameters -np ${no_mpi_proc} "./utest-Poisson")
add_executable(utest-Poisson utest-Poisson.cpp)
add_test(NAME utest-Poisson COMMAND "mpirun" ${utest-Poisson_parameters} )
target_link_libraries(utest-Poisson common mesh interpolation linear_system boost_unit_test_framework)

set(utest-weak-dir-bc_parameters -np ${no_mpi_proc} "./utest-weak-dir-bc")
add_executable(utest-weak-dir-bc utest-weak-dir-bc.cpp)
add_test(NAME utest-weak-dir-bc COMMAND "mpirun" ${utest-weak-dir-bc_parameters} )
target_link_libraries(utest-weak-dir-bc common mesh interpolation linear_system fe boost_unit_test_framework)


# Copy the meshes from the mesh repository into the build tree
add_custom_command(TARGET utest-Poisson
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_2D_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_2D_tri_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_2D_tri_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/weak_bc_square_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/weak_bc_square_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/weak_bc_square_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/weak_bc_square_p4.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_3D_tet_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_3D_tet_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_3D_tet_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Poisson/L_shape_3D_tet_p4.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR})

endif( PDEKIT_HAVE_TRILINOS )
