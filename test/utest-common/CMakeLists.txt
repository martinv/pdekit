# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

set( no_mpi_proc 2 )

add_executable(utest-sizes-common utest-sizes-common.cpp)
add_test(utest-sizes-common utest-sizes-common)

add_executable(utest-core utest-core.cpp)
add_test(utest-core utest-core)
target_link_libraries(utest-core common boost_unit_test_framework)

add_executable(utest-array-view utest-array-view.cpp)
add_test(utest-array-view utest-array-view)
target_link_libraries(utest-array-view boost_unit_test_framework)

add_executable(utest-blockarray utest-blockarray.cpp)
add_test(utest-blockarray utest-blockarray)
target_link_libraries(utest-blockarray boost_unit_test_framework)

add_executable(utest-blockmultiarray utest-blockmultiarray.cpp)
add_test(utest-blockmultiarray utest-blockmultiarray)
target_link_libraries(utest-blockmultiarray boost_unit_test_framework)

add_executable(utest-component utest-component.cpp)
add_test(utest-component utest-component)
target_link_libraries(utest-component common boost_unit_test_framework)

add_executable(utest-hash utest-hash.cpp)
add_test(utest-hash utest-hash)
target_link_libraries(utest-hash common boost_unit_test_framework)

add_executable(utest-option-map utest-option-map.cpp)
add_test(utest-option-map utest-option-map)
target_link_libraries(utest-option-map common boost_unit_test_framework)

add_executable(utest-range utest-range.cpp)
add_test(utest-range utest-range)
target_link_libraries(utest-range common boost_unit_test_framework)

add_executable(utest-singleton utest-singleton.cpp)
add_test(utest-singleton utest-singleton)
target_link_libraries(utest-singleton common boost_unit_test_framework)

add_executable(utest-factory utest-factory.cpp)
add_test(utest-factory utest-factory)
target_link_libraries(utest-factory common boost_unit_test_framework)

add_executable(utest-tags utest-tags.cpp)
add_test(utest-tags utest-tags)
target_link_libraries(utest-tags common boost_unit_test_framework)

add_executable(utest-uri utest-uri.cpp)
add_test(utest-uri utest-uri)
target_link_libraries(utest-uri common boost_unit_test_framework)

add_executable(utest-vector-utils utest-vector-utils.cpp)
add_test(utest-vector-utils utest-vector-utils)
target_link_libraries(utest-vector-utils common boost_unit_test_framework)

add_executable(check_mpi_functionality check_mpi_functionality.cpp)
target_link_libraries(check_mpi_functionality common boost_unit_test_framework)

add_executable(utest-mpi-send-recv utest-mpi-send-recv.cpp)
add_test(NAME utest-mpi-send-recv COMMAND "mpiexec" -np ${no_mpi_proc} "./utest-mpi-send-recv")
target_link_libraries(utest-mpi-send-recv common boost_unit_test_framework)
