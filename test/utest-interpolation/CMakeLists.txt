
add_executable(utest-fe-values utest-fe-values.cpp)
add_test(utest-fe-values utest-fe-values)
target_link_libraries(utest-fe-values mesh interpolation boost_unit_test_framework)

add_executable(utest-cell-distance utest-cell-distance.cpp)
add_test(utest-cell-distance utest-cell-distance)
target_link_libraries(utest-cell-distance interpolation boost_unit_test_framework)

add_executable(utest-function-space utest-function-space.cpp)
add_test(utest-function-space utest-function-space)
target_link_libraries(utest-function-space mesh interpolation boost_unit_test_framework)

add_executable(utest-element-adapt-interpolator utest-element-adapt-interpolator.cpp)
add_test(utest-element-adapt-interpolator utest-element-adapt-interpolator)
target_link_libraries(utest-element-adapt-interpolator mesh interpolation boost_unit_test_framework)

add_executable(utest-geometry-cache utest-geometry-cache.cpp)
add_test(utest-geometry-cache utest-geometry-cache)
target_link_libraries(utest-geometry-cache common mesh interpolation boost_unit_test_framework)

add_executable(utest-geometry-metric utest-geometry-metric.cpp)
add_test(utest-geometry-metric utest-geometry-metric)
target_link_libraries(utest-geometry-metric common mesh interpolation boost_unit_test_framework)

add_executable(utest-solution-space-metric utest-solution-space-metric.cpp)
add_test(utest-solution-space-metric utest-solution-space-metric)
target_link_libraries(utest-solution-space-metric common mesh interpolation boost_unit_test_framework)

add_executable(utest-polynomial-metric utest-polynomial-metric.cpp)
add_test(utest-polynomial-metric utest-polynomial-metric)
target_link_libraries(utest-polynomial-metric common mesh interpolation boost_unit_test_framework)

add_executable(utest-flux-space-metric utest-flux-space-metric.cpp)
add_test(utest-flux-space-metric utest-flux-space-metric)
target_link_libraries(utest-flux-space-metric common mesh interpolation physics boost_unit_test_framework)

add_executable(utest-operator-residual utest-operator-residual.cpp)
add_test(utest-operator-residual utest-operator-residual)
target_link_libraries(utest-operator-residual common mesh interpolation physics boost_unit_test_framework)

add_executable(utest-scalar-mesh-function utest-scalar-mesh-function.cpp)
add_test(utest-scalar-mesh-function utest-scalar-mesh-function)
target_link_libraries(utest-scalar-mesh-function common mesh interpolation boost_unit_test_framework)

add_executable(utest-vector-mesh-function utest-vector-mesh-function.cpp)
add_test(utest-vector-mesh-function utest-vector-mesh-function)
target_link_libraries(utest-vector-mesh-function common mesh interpolation boost_unit_test_framework)

add_executable(utest-mesh-function-ops utest-mesh-function-ops.cpp)
add_test(utest-mesh-function-ops utest-mesh-function-ops)
target_link_libraries(utest-mesh-function-ops common mesh interpolation boost_unit_test_framework)

add_custom_command(TARGET utest-geometry-metric
                POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p1_regular.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p2_regular.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mixed_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mixed_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_tri_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_cube_tet_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_cube_tet_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_cube_tet_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_circle_tri_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_circle_tri_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_circle_tri_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_sphere_tet_p1.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_sphere_tet_p2.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_disc/unit_sphere_tet_p3.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
               )
