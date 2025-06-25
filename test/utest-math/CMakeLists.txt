# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

add_executable(utest-jacobi-polynomial utest-jacobi-polynomial.cpp)
add_test(utest-jacobi-polynomial utest-jacobi-polynomial)
target_link_libraries(utest-jacobi-polynomial math boost_unit_test_framework)

add_executable(utest-et-dynamic-matrix-vector-ops utest-et-dynamic-matrix-vector-ops.cpp)
add_test(utest-et-dynamic-matrix-vector-ops utest-et-dynamic-matrix-vector-ops)
target_link_libraries(utest-et-dynamic-matrix-vector-ops math)

add_executable(utest-et-static-matrix-vector-ops utest-et-static-matrix-vector-ops.cpp)
add_test(utest-et-static-matrix-vector-ops utest-et-static-matrix-vector-ops)
target_link_libraries(utest-et-static-matrix-vector-ops math)

add_executable(utest-et-corectness utest-et-corectness.cpp)
add_test(utest-et-corectness utest-et-corectness)
target_link_libraries(utest-et-corectness math boost_unit_test_framework)

add_executable(utest-sparse-matrix utest-sparse-matrix.cpp)
add_test(utest-sparse-matrix utest-sparse-matrix)
target_link_libraries(utest-sparse-matrix math boost_unit_test_framework)

add_executable(utest-block-sparse-matrix utest-block-sparse-matrix.cpp)
add_test(utest-block-sparse-matrix utest-block-sparse-matrix)
target_link_libraries(utest-block-sparse-matrix math boost_unit_test_framework)

add_executable(utest-matrix-inverter utest-matrix-inverter.cpp)
add_test(utest-matrix-inverter utest-matrix-inverter)
target_link_libraries(utest-matrix-inverter math boost_unit_test_framework)

add_executable(utest-matrix-decompositions utest-matrix-decompositions.cpp)
add_test(utest-matrix-decompositions utest-matrix-decompositions)
target_link_libraries(utest-matrix-decompositions math boost_unit_test_framework)

add_executable(utest-mat-array utest-mat-array.cpp)
add_test(utest-mat-array utest-mat-array)
target_link_libraries(utest-mat-array math boost_unit_test_framework)

add_executable(utest-matrix-condition-number utest-matrix-condition-number.cpp)
add_test(utest-matrix-condition-number utest-matrix-condition-number)
target_link_libraries(utest-matrix-condition-number math boost_unit_test_framework)

add_executable(utest-polylib utest-polylib.cpp)
add_test(utest-polylib utest-polylib)
target_link_libraries(utest-polylib math boost_unit_test_framework)

add_executable(utest-polynomial utest-polynomial.cpp)
add_test(utest-polynomial utest-polynomial)
target_link_libraries(utest-polynomial math boost_unit_test_framework)

add_executable(utest-gram-schmidt utest-gram-schmidt.cpp)
add_test(utest-gram-schmidt utest-gram-schmidt)
target_link_libraries(utest-gram-schmidt math boost_unit_test_framework)
