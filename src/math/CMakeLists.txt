#CmakeLists.txt in math dir

include(PDEKITFiles.cmake)
include(traits/PDEKITFiles.cmake)
include(polynomials/PDEKITFiles.cmake)
include(unary_ops/PDEKITFiles.cmake )
include(binary_ops/PDEKITFiles.cmake)
include(decompositions/PDEKITFiles.cmake)
include(algo/PDEKITFiles.cmake)

add_library( math SHARED ${PDEKIT_Math_HEADERS} 
                         ${PDEKIT_Math_Traits_HEADERS}
                         ${PDEKIT_Math_Binary_Ops_HEADERS}
                         ${PDEKIT_Math_Polynomials_HEADERS}
                         ${PDEKIT_Math_Polynomials_SOURCES}
                         ${PDEKIT_Math_Vector_Unary_Ops_HEADERS}
                         ${PDEKIT_Math_Matrix_Unary_Ops_HEADERS}
                         ${PDEKIT_Math_Matrix_Decompositions_HEADERS}
                         ${PDEKIT_Math_Algo_HEADERS})

if ( PDEKIT_USE_INTERNAL_BLAS_LAPACK )
  add_dependencies(math InternalLapack)
endif ( PDEKIT_USE_INTERNAL_BLAS_LAPACK )

target_compile_options(math PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(math ${BLAS_LAPACK_LIBRARIES})

install(FILES ${PDEKIT_Math_HEADERS} DESTINATION include/math)
install(FILES ${PDEKIT_Math_Traits_HEADERS} DESTINATION include/math/traits)
install(FILES ${PDEKIT_Math_Binary_Ops_HEADERS} DESTINATION include/math/binary_ops)
install(TARGETS math DESTINATION lib)

