# Include a list of files in the directory io/gmsh
include(io/gmsh/PDEKITFiles.cmake)

# Include a list of files in the directory io/vtk
include(io/vtk/PDEKITFiles.cmake)

# List of files directly present in the directory mesh/io
list(APPEND PDEKIT_Mesh_Generic_IO_HEADERS
 io/MeshIO.hpp
 io/MeshCreator.hpp
 io/MeshManipulator.hpp
)

# List of files directly present in the directory mesh/io
list(APPEND PDEKIT_Mesh_Generic_IO_SOURCES
 io/MeshCreator.cpp
 io/MeshManipulator.cpp
)

list(APPEND PDEKIT_Mesh_IO_HEADERS
 ${PDEKIT_Mesh_Generic_IO_HEADERS}
 ${PDEKIT_Mesh_Gmsh_HEADERS} 
 ${PDEKIT_Mesh_VTK_HEADERS} 
)

list(APPEND PDEKIT_Mesh_IO_SOURCES
 ${PDEKIT_Mesh_Generic_IO_SOURCES}
 ${PDEKIT_Mesh_Gmsh_SOURCES} 
 ${PDEKIT_Mesh_CF_SOURCES} 
 ${PDEKIT_Mesh_VTK_SOURCES}
)

