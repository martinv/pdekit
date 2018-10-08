# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

#add_executable(utest-mesh-entity utest-mesh-entity.cpp)
#add_test(utest-mesh-entity utest-mesh-entity)
#target_link_libraries(utest-mesh-entity common mesh boost_unit_test_framework)

add_executable(utest-quadrature-type utest-quadrature-type.cpp)
add_test(utest-quadrature-type utest-quadrature-type)
target_link_libraries(utest-quadrature-type common mesh boost_unit_test_framework)

add_executable(utest-quadrature-permutation utest-quadrature-permutation.cpp)
add_test(utest-quadrature-permutation utest-quadrature-permutation)
target_link_libraries(utest-quadrature-permutation common mesh boost_unit_test_framework)

add_executable(utest-quadrature-point-set utest-quadrature-point-set.cpp)
add_test(utest-quadrature-point-set utest-quadrature-point-set)
target_link_libraries(utest-quadrature-point-set common math mesh boost_unit_test_framework)

add_executable(utest-quadrature-factory utest-quadrature-factory.cpp)
add_test(utest-quadrature-factory utest-quadrature-factory)
target_link_libraries(utest-quadrature-factory common mesh boost_unit_test_framework)

add_executable(utest-quadrature-adapt-transform-algo utest-quadrature-adapt-transform-algo.cpp)
add_test(utest-quadrature-adapt-transform-algo utest-quadrature-adapt-transform-algo)
target_link_libraries(utest-quadrature-adapt-transform-algo common mesh boost_unit_test_framework)
