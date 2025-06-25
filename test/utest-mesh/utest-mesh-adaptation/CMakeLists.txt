add_executable(utest-adapt-mesh utest-adapt-mesh.cpp)
add_test(utest-adapt-mesh utest-adapt-mesh)
target_link_libraries(utest-adapt-mesh common mesh boost_unit_test_framework)

add_executable(utest-local-interpolator utest-local-interpolator.cpp)
add_test(utest-local-interpolator utest-local-interpolator)
target_link_libraries(utest-local-interpolator common mesh boost_unit_test_framework)

add_executable(utest-cell-split-strategy utest-cell-split-strategy.cpp)
add_test(utest-cell-split-strategy utest-cell-split-strategy)
target_link_libraries(utest-cell-split-strategy common mesh boost_unit_test_framework)

add_custom_command(TARGET utest-adapt-mesh
                  POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/test_p1_tri_no_bdry.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                  POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PDEKIT_MESH_REPO_DIR}/unit_interval/unit_square_mini.msh ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}
                  )



