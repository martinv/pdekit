#CmakeLists.txt in Mesh dir

#Adds a library called Mesh (libmesh.a under Linux) from the source files

include(PDEKITFiles.cmake)
include(containers/PDEKITFiles.cmake)
include(iterators/PDEKITFiles.cmake)
include(algorithm/PDEKITFiles.cmake)
include(view/PDEKITFiles.cmake)
include(std_region/PDEKITFiles.cmake)
include(local_topology/PDEKITFiles.cmake)
include(point_set/PDEKITFiles.cmake)
include(shape_function/PDEKITFiles.cmake)
include(io/PDEKITFiles.cmake)
include(adaptation/PDEKITFiles.cmake)

add_library(mesh SHARED ${PDEKIT_Mesh_HEADERS} ${PDEKIT_Mesh_SOURCES} 
                        ${PDEKIT_Mesh_Containers_HEADERS} ${PDEKIT_Mesh_Containers_SOURCES}
                        ${PDEKIT_Mesh_Iterators_HEADERS} ${PDEKIT_Mesh_Iterators_SOURCES}
                        ${PDEKIT_Mesh_Algorithm_HEADERS} ${PDEKIT_Mesh_Algorithm_SOURCES}
                        ${PDEKIT_Mesh_View_HEADERS} ${PDEKIT_Mesh_View_SOURCES}
                        ${PDEKIT_Mesh_IO_HEADERS} ${PDEKIT_Mesh_IO_SOURCES} 
                        ${PDEKIT_Mesh_Std_Region_HEADERS} ${PDEKIT_Mesh_Std_Region_SOURCES} 
                        ${PDEKIT_Mesh_Local_Topology_HEADERS} ${PDEKIT_Mesh_Local_Topology_SOURCES} 
                        ${PDEKIT_Mesh_Shape_Function_HEADERS} ${PDEKIT_Mesh_Shape_Function_SOURCES}
                        ${PDEKIT_Mesh_Point_Set_HEADERS} ${PDEKIT_Mesh_Point_Set_SOURCES}
                        ${PDEKIT_Mesh_Adaptation_HEADERS} ${PDEKIT_Mesh_Adaptation_SOURCES})

target_compile_options(mesh PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})
target_link_libraries(mesh common math)

install(FILES ${PDEKIT_Mesh_HEADERS} DESTINATION include/mesh)
install(FILES ${PDEKIT_Mesh_Containers_HEADERS} DESTINATION include/mesh/containers)
install(FILES ${PDEKIT_Mesh_Iterators_HEADERS} DESTINATION include/mesh/iterators)
install(FILES ${PDEKIT_Mesh_Algorithm_HEADERS} DESTINATION include/mesh/algorithm)
install(FILES ${PDEKIT_Mesh_View_HEADERS} DESTINATION include/mesh/view)
install(FILES ${PDEKIT_Mesh_IO_HEADERS} DESTINATION include/mesh/io)
install(FILES ${PDEKIT_Mesh_Std_Region_HEADERS} DESTINATION include/mesh/std_region)
install(FILES ${PDEKIT_Mesh_Local_Topology_HEADERS} DESTINATION include/mesh/local_topology)
install(FILES ${PDEKIT_Mesh_Shape_Function_HEADERS} DESTINATION include/mesh/shape_function)
install(FILES ${PDEKIT_Mesh_Point_Set_HEADERS} DESTINATION include/mesh/point_set)
install(FILES ${PDEKIT_Mesh_Adaptation_HEADERS} DESTINATION include/mesh/adaptation)

install(TARGETS mesh DESTINATION lib)
