# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

add_executable(utest-graph utest-graph.cpp)
add_test(utest-graph utest-graph)
target_link_libraries(utest-graph common boost_unit_test_framework)

add_executable(utest-graph-partitioner utest-graph-partitioner.cpp)
add_test(utest-graph-partitioner utest-graph-partitioner)
target_link_libraries(utest-graph-partitioner common graph boost_unit_test_framework)
