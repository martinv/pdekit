# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

if ( PDEKIT_HAVE_TRILINOS )

add_executable(utest-helmholtz-hdg utest-helmholtz-hdg.cpp)
add_test(utest-helmholtz-hdg utest-helmholtz-hdg)
target_link_libraries(utest-helmholtz-hdg common graph mesh interpolation solver linear_system boost_unit_test_framework)

add_executable(utest-helmholtz-cg-hdg utest-helmholtz-cg-hdg.cpp)
add_test(utest-helmholtz-cg-hdg utest-helmholtz-cg-hdg)
target_link_libraries(utest-helmholtz-cg-hdg common graph mesh interpolation solver linear_system boost_unit_test_framework)

# Copy the meshes from the mesh repository into the build tree
add_custom_command(TARGET utest-helmholtz-cg-hdg
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
