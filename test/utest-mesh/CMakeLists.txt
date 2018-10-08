# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

add_executable(utest-cell-buffer utest-cell-buffer.cpp)
add_test(utest-cell-buffer utest-cell-buffer)
target_link_libraries(utest-cell-buffer common mesh boost_unit_test_framework)

add_executable(utest-cell-path utest-cell-path.cpp)
add_test(utest-cell-path utest-cell-path)
target_link_libraries(utest-cell-path common mesh boost_unit_test_framework)

add_executable(utest-point utest-point.cpp)
add_test(utest-point utest-point)
target_link_libraries(utest-point common mesh boost_unit_test_framework)

add_executable(utest-cell-adapt-op utest-cell-adapt-op.cpp)
add_test(utest-cell-adapt-op utest-cell-adapt-op)
target_link_libraries(utest-cell-adapt-op common mesh boost_unit_test_framework)

add_executable(utest-cell-geometry utest-cell-geometry.cpp)
add_test(utest-cell-geometry utest-cell-geometry)
target_link_libraries(utest-cell-geometry common mesh boost_unit_test_framework)

add_executable(utest-entity-hash-collision utest-entity-hash-collision.cpp)
add_test(utest-entity-hash-collision utest-entity-hash-collision)
target_link_libraries(utest-entity-hash-collision common mesh boost_unit_test_framework)

add_executable(utest-sizes-mesh utest-sizes-mesh.cpp)
add_test(utest-sizes-mesh utest-sizes-mesh)
target_link_libraries(utest-sizes-mesh common mesh)

add_executable(utest-mesh-geometry utest-mesh-geometry.cpp)
add_test(utest-mesh-geometry utest-mesh-geometry)
target_link_libraries(utest-mesh-geometry common mesh boost_unit_test_framework)

add_executable(utest-element-type utest-element-type.cpp)
add_test(utest-element-type utest-element-type)
target_link_libraries(utest-element-type common mesh)

add_executable(utest-reference-topology utest-reference-topology.cpp)
add_test(utest-reference-topology utest-reference-topology)
target_link_libraries(utest-reference-topology common mesh boost_unit_test_framework)

add_executable(utest-dof-point-set utest-dof-point-set.cpp)
add_test(utest-dof-point-set utest-dof-point-set)
target_link_libraries(utest-dof-point-set common mesh boost_unit_test_framework)

add_executable(utest-std-region-data-map utest-std-region-data-map.cpp)
add_test(utest-std-region-data-map utest-std-region-data-map)
target_link_libraries(utest-std-region-data-map common mesh boost_unit_test_framework)

add_executable(utest-mesh-entity utest-mesh-entity.cpp)
add_test(utest-mesh-entity utest-mesh-entity)
target_link_libraries(utest-mesh-entity common mesh boost_unit_test_framework)

add_executable(utest-mesh-manipulator utest-mesh-manipulator.cpp)
add_test(utest-mesh-manipulator utest-mesh-manipulator)
target_link_libraries(utest-mesh-manipulator common mesh interpolation linear_system boost_unit_test_framework)

add_executable(utest-cell-subdomain-tuple utest-cell-subdomain-tuple.cpp)
add_test(utest-cell-subdomain-tuple utest-cell-subdomain-tuple)
target_link_libraries(utest-cell-subdomain-tuple common mesh boost_unit_test_framework)

add_executable(utest-topology-containers utest-topology-containers.cpp)
add_test(utest-topology-containers utest-topology-containers)
target_link_libraries(utest-topology-containers common mesh boost_unit_test_framework)

add_executable(utest-cell-topology-iterator utest-cell-topology-iterator.cpp)
add_test(utest-cell-topology-iterator utest-cell-topology-iterator)
target_link_libraries(utest-cell-topology-iterator mesh boost_unit_test_framework)

add_executable(utest-topology-bdry-iterator utest-topology-bdry-iterator.cpp)
add_test(utest-topology-bdry-iterator utest-topology-bdry-iterator)
target_link_libraries(utest-topology-bdry-iterator mesh boost_unit_test_framework)

add_executable(utest-skeleton-iterator utest-skeleton-iterator.cpp)
add_test(utest-skeleton-iterator utest-skeleton-iterator)
target_link_libraries(utest-skeleton-iterator mesh boost_unit_test_framework)

add_executable(utest-dof-iterator utest-dof-iterator.cpp)
add_test(utest-dof-iterator utest-dof-iterator)
target_link_libraries(utest-dof-iterator mesh boost_unit_test_framework)

add_executable(utest-partition-mesh utest-partition-mesh.cpp)
add_test(utest-partition-mesh utest-partition-mesh)
target_link_libraries(utest-partition-mesh common graph mesh boost_unit_test_framework)

add_executable(utest-gmsh-reader-writer utest-gmsh-reader-writer.cpp)
add_test(utest-gmsh-reader-writer utest-gmsh-reader-writer)
target_link_libraries(utest-gmsh-reader-writer common mesh boost_unit_test_framework)

add_executable(utest-vtk-reader-writer utest-vtk-reader-writer.cpp)
add_test(utest-vtk-reader-writer utest-vtk-reader-writer)
target_link_libraries(utest-vtk-reader-writer common mesh boost_unit_test_framework)

add_executable(utest-reference-element utest-reference-element.cpp)
add_test(utest-reference-element utest-reference-element)
target_link_libraries(utest-reference-element common mesh boost_unit_test_framework)

add_executable(utest-topology-algorithms utest-topology-algorithms.cpp)
add_test(utest-topology-algorithms utest-topology-algorithms)
target_link_libraries(utest-topology-algorithms common mesh boost_unit_test_framework)

add_executable(utest-topology-predicates utest-topology-predicates.cpp)
add_test(utest-topology-predicates utest-topology-predicates)
target_link_libraries(utest-topology-predicates common mesh boost_unit_test_framework)

add_executable(utest-upgrade-mesh utest-upgrade-mesh.cpp)
add_test(utest-upgrade-mesh utest-upgrade-mesh)
target_link_libraries(utest-upgrade-mesh common mesh boost_unit_test_framework)

add_executable(utest-geo-transform-mesh utest-geo-transform-mesh.cpp)
add_test(utest-geo-transform-mesh utest-geo-transform-mesh)
target_link_libraries(utest-geo-transform-mesh common mesh boost_unit_test_framework)

add_executable(utest-mesh-creator utest-mesh-creator.cpp)
add_test(utest-mesh-creator utest-mesh-creator)
target_link_libraries(utest-mesh-creator common mesh boost_unit_test_framework)

# =============================================================================

add_subdirectory(shape_function)
add_subdirectory(utest-mesh-adaptation)

# =============================================================================


# Copy the meshes from the mesh repository into the build tree
add_custom_command(TARGET utest-mesh-geometry
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/cube.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )

add_custom_command(TARGET utest-gmsh-reader-writer
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p1_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p2_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p3_tri.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/sinus_bump_3d/bump_p1_tet.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )

add_custom_command(TARGET utest-vtk-reader-writer
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/Cube.vtk ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/Leaflet_polygon.vtk ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                )

add_custom_command(TARGET utest-dof-iterator
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/rectangle_mixed_elem_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/rectangle_mixed_elem_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/rectangle_mixed_elem_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_cube_mixed_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )

add_custom_command(TARGET utest-skeleton-iterator
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mini.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )

add_custom_command(TARGET utest-topology-algorithms
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p1_tri_no_bdry.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )

#COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CF_RESOURCE_DIR}/rotation-qd-p1.neu ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
# ADD_EXECUTABLE(SFTest SFTest.cpp)
# ADD_TEST(ShapeFunctionTest SFTest) # ADD_TEST(TestName ExeName)


