
add_executable(linear-advection-continuous-rds-2D-implicit LinearAdvection-Continuous-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(linear-advection-continuous-rds-2D-implicit common math mesh rdm)

add_executable(advection-continuous-rds-2D-explicit Advection-Continuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-continuous-rds-2D-explicit common math mesh rdm time)

add_executable(advection-rds-2D Advection-RDS-2D.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-rds-2D common math mesh rdm time)

add_executable(advection-rds-2D-spectral Advection-RDS-2D-Spectral.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-rds-2D-spectral common math mesh rdm time)

add_executable(advection-continuous-rds-2D-implicit Advection-Continuous-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-continuous-rds-2D-implicit common math mesh rdm)

add_executable(diffusion-continuous-rds-2D-explicit Diffusion-Continuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(diffusion-continuous-rds-2D-explicit common math mesh rdm time)

add_executable(advection-src-term-continuous-rds-2D Advection-SrcTerm-Continuous-RDS-2D.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-src-term-continuous-rds-2D common math mesh rdm)

add_executable(advection-multigrid-continuous-rds-2D Advection-Multigrid-Continuous-RDS-2D.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-multigrid-continuous-rds-2D common math mesh solver rdm time)

add_executable(advection-discontinuous-rds-2D-explicit Advection-Discontinuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-discontinuous-rds-2D-explicit common math mesh rdm)

add_executable(advection-discontinuous-rds-2D-implicit Advection-Discontinuous-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-discontinuous-rds-2D-implicit common math mesh rdm)

add_executable(advection-rds-2D-explicit-h-adapt Advection-RDS-2D-Explicit-h-Adapt.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-rds-2D-explicit-h-adapt common math mesh rdm)

add_executable(advection-rds-2D-explicit-p-adapt Advection-RDS-2D-Explicit-p-Adapt.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-rds-2D-explicit-p-adapt common math mesh rdm)

add_executable(advection-continuous-rds-3D-explicit Advection-Continuous-RDS-3D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-continuous-rds-3D-explicit common math mesh rdm time)

add_executable(advection-rds-3D-implicit Advection-RDS-3D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(advection-rds-3D-implicit common math mesh rdm)

add_executable(Burgers-rds-2D Burgers-RDS-2D.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(Burgers-rds-2D common math mesh rdm)

add_executable(Burgers-rds-2D-adapt Burgers-RDS-2D-Adapt.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(Burgers-rds-2D-adapt common math mesh rdm)

add_executable(jets-continuous-rds-2D-explicit Jets-Continuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(jets-continuous-rds-2D-explicit common math mesh rdm)

add_executable(jets-continuous-rds-2D-implicit Jets-Continuous-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(jets-continuous-rds-2D-implicit common math mesh rdm)

add_executable(jets-discontinuous-rds-2D-explicit Jets-Discontinuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(jets-discontinuous-rds-2D-explicit common math mesh rdm)

add_executable(jets3d-rds Jets3D-RDS.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(jets3d-rds common math mesh rdm)

add_executable(bump-rds-2D-explicit Bump-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(bump-rds-2D-explicit common math mesh rdm)

add_executable(bump-rds-2D-implicit Bump-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(bump-rds-2D-implicit common math mesh rdm)

add_executable(bump-rds-3D-explicit Bump-RDS-3D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(bump-rds-3D-explicit common math mesh rdm)

add_executable(bump-rds-3D-implicit Bump-RDS-3D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(bump-rds-3D-implicit common math mesh rdm)

add_executable(ringleb-2D Ringleb.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(ringleb-2D common math mesh rdm)

add_executable(naca-continuous-rds-2d-explicit Naca-Continuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(naca-continuous-rds-2d-explicit common math mesh rdm time)

add_executable(naca-discontinuous-rds-2d-explicit Naca-Discontinuous-RDS-2D-Explicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(naca-discontinuous-rds-2d-explicit common math mesh rdm)

add_executable(naca-rds-2d-implicit Naca-RDS-2D-Implicit.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(naca-rds-2d-implicit common math mesh rdm)

add_executable(naca-rds-2d-adapt Naca-RDS-2D-Adapt.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(naca-rds-2d-adapt common math mesh rdm)

add_executable(naca2D-multigrid-continuous-rds Naca2D-Multigrid-Continuous-RDS.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(naca2D-multigrid-continuous-rds common math mesh solver rdm time)

add_executable(btc0-rds BTC0-RDS.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(btc0-rds common math mesh rdm)

add_executable(oneraM6-rds OneraM6-RDS.cpp PredefinedFunctions.hpp PredefinedFunctions.cpp)
target_link_libraries(oneraM6-rds common math mesh rdm)

add_executable(supersonic-cylinder-2D SupersonicCylinder-2D.cpp)
target_link_libraries(supersonic-cylinder-2D common math mesh rdm)

add_executable(supersonic-cylinder-2D-adapt SupersonicCylinder-2D-Adapt.cpp)
target_link_libraries(supersonic-cylinder-2D-adapt common math mesh rdm)

add_executable(subsonic-cylinder-rds-2D SubsonicCylinder-RDS-2D.cpp)
target_link_libraries(subsonic-cylinder-rds-2D common math mesh rdm)

add_executable(postprocess-data PostprocessData.cpp)
target_link_libraries(postprocess-data common math mesh rdm)

# Copy the meshes from the mesh repository into the build tree
add_custom_command(TARGET linear-advection-continuous-rds-2D-implicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET advection-continuous-rds-2D-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/advection_2d/advection-tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/advection_2d/advection-quad.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/advection_2d/advection-delquad.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET diffusion-continuous-rds-2D-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET Burgers-rds-2D
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Burgers/square_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Burgers/square_tri_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Burgers/square_tri_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Burgers/square_tri_p4.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET advection-continuous-rds-3D-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/advection_3d/advection-tet.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/advection_3d/advection-hex.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET jets-continuous-rds-2D-explicit
		POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/trapezoid/trapezoid-tri-p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFD_INTDIR}
)

add_custom_command(TARGET jets3d-rds
		POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/trapezoid/trapezoid-tet-p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFD_INTDIR}
)

add_custom_command(TARGET bump-rds-2D-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/rectangle_channel_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_2d/bump_p1_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_2d/bump_p2_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET bump-rds-3D-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_3d/bump_p1_tet.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_3d/bump_p2_tet.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_3d/bump_p3_tet.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET ringleb-2D
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/ringleb/ringleb_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR})

add_custom_command(TARGET naca-continuous-rds-2d-explicit
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/naca0012_p1_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/naca0012_p2_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/naca0012_p3_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/mixed_elem_mesh/naca0012_p1_mixed_elem.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/mixed_elem_mesh/naca0012_p2_mixed_elem.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/mixed_elem_mesh/naca0012_p3_mixed_elem.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/naca0012_2d/mixed_elem_mesh/naca0012_p4_mixed_elem.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}

)

add_custom_command(TARGET btc0-rds
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/cigar/btc0-mesh-p2-coarse.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET oneraM6-rds
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/onera_m6/M6_unstr_euler_c1_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/onera_m6/M6_unstr_euler_c2_p2_rescaled.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET supersonic-cylinder-2D
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/supersonic_cylinder/supersonic_cylinder_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/supersonic_cylinder/supersonic_cylinder_tri_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/supersonic_cylinder/supersonic_cylinder_tri_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET subsonic-cylinder-rds-2D
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/subsonic_cylinder/cylinder2_regular_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/subsonic_cylinder/cylinder2_regular_p2p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/subsonic_cylinder/cylinder2_regular_p3p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET advection-src-term-continuous-rds-2D
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)

add_custom_command(TARGET postprocess-data
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/onera_m6/M6_unstr_euler_c2_p2_sample_result.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
)
