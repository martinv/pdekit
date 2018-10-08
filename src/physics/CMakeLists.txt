#CmakeLists.txt in physics dir

#Adds a library called physics (libphysics.{a,so} under Linux) from the source files

include(PDEKITFiles.cmake)
include(scalar/PDEKITFiles.cmake)
include(euler/PDEKITFiles.cmake)

add_library(physics SHARED ${PDEKIT_Physics_HEADERS} 
                           ${PDEKIT_Physics_Scalar_Advection_HEADERS} ${PDEKIT_Physics_Scalar_Advection_SOURCES} 
                           ${PDEKIT_Physics_Euler_HEADERS} ${PDEKIT_Physics_Euler_SOURCES})

target_compile_options(physics PUBLIC ${PDEKIT_CMAKE_CXX_FLAGS})

install(FILES ${PDEKIT_Physics_HEADERS} DESTINATION include/physics)
install(FILES ${PDEKIT_Physics_Scalar_Advection_HEADERS} DESTINATION include/physics/scalar)
install(FILES ${PDEKIT_Physics_Euler_HEADERS} DESTINATION include/physics/euler)


