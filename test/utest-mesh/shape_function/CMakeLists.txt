# ADD_TEST(TestName ExeName)
# ADD_CUSTOM_TARGET(check COMMAND ${CMAKE_CTEST_COMMAND}
#                   DEPENDS ExeName)

add_executable(utest-LagrangeSF utest-LagrangeSF.cpp
               fixture/ReferenceShapeFunctionsTetra.hpp
               fixture/ReferenceShapeFunctionsTetra.cpp)
add_test(utest-LagrangeSF utest-LagrangeSF)
target_link_libraries(utest-LagrangeSF common math mesh interpolation boost_unit_test_framework)

add_executable(utest-Dubiner-expansion-simplex utest-Dubiner-expansion-simplex.cpp 
					       fixture/ReferenceShapeFunctionsTriag.hpp
					       fixture/ReferenceShapeFunctionsTriag.cpp
					       fixture/ReferenceShapeFunctionsTetra.hpp
					       fixture/ReferenceShapeFunctionsTetra.cpp)
add_test(utest-Dubiner-expansion-simplex utest-Dubiner-expansion-simplex)
target_link_libraries(utest-Dubiner-expansion-simplex common mesh boost_unit_test_framework)

add_executable(utest-Carnevali-expansion-simplex utest-Carnevali-expansion-simplex.cpp 
					       fixture/ReferenceShapeFunctionsTriag.hpp
					       fixture/ReferenceShapeFunctionsTriag.cpp
					       fixture/ReferenceShapeFunctionsTetra.hpp
					       fixture/ReferenceShapeFunctionsTetra.cpp)
add_test(utest-Carnevali-expansion-simplex utest-Carnevali-expansion-simplex)
target_link_libraries(utest-Carnevali-expansion-simplex common mesh boost_unit_test_framework)

add_executable(utest-modal-expansion-line utest-modal-expansion-line.cpp 
					       fixture/ReferenceShapeFunctionsLine.hpp
					       fixture/ReferenceShapeFunctionsLine.cpp)
add_test(utest-modal-expansion-line utest-modal-expansion-line)
target_link_libraries(utest-modal-expansion-line common mesh boost_unit_test_framework)

add_executable(utest-modal-expansion-quad-hexa utest-modal-expansion-quad-hexa.cpp 
					       fixture/ReferenceShapeFunctionsQuad.hpp
                 fixture/ReferenceShapeFunctionsQuad.cpp
                 fixture/ReferenceShapeFunctionsHexa.hpp
                 fixture/ReferenceShapeFunctionsHexa.cpp)
add_test(utest-modal-expansion-quad-hexa utest-modal-expansion-quad-hexa)
target_link_libraries(utest-modal-expansion-quad-hexa common mesh boost_unit_test_framework)

add_executable(utest-ortho-basis-gram-schmidt utest-ortho-basis-gram-schmidt.cpp)
add_test(utest-ortho-basis-gram-schmidt utest-ortho-basis-gram-schmidt)
target_link_libraries(utest-ortho-basis-gram-schmidt common mesh boost_unit_test_framework)
