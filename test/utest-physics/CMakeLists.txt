
add_executable(utest-linearadv2d-model utest-linearadv2d-model.cpp)
add_test(utest-linearadv2d-model utest-linearadv2d-model)
target_link_libraries(utest-linearadv2d-model physics boost_unit_test_framework)

add_executable(utest-adv-diff-fos-model-2D utest-adv-diff-fos-model-2D.cpp)
add_test(utest-adv-diff-fos-model-2D utest-adv-diff-fos-model-2D)
target_link_libraries(utest-adv-diff-fos-model-2D physics boost_unit_test_framework)

add_executable(utest-euler-model-2D utest-euler-model-2D.cpp)
add_test(utest-euler-model-2D utest-euler-model-2D)
target_link_libraries(utest-euler-model-2D physics boost_unit_test_framework)

add_executable(utest-euler-model-3D utest-euler-model-3D.cpp)
add_test(utest-euler-model-3D utest-euler-model-3D)
target_link_libraries(utest-euler-model-3D physics boost_unit_test_framework)
