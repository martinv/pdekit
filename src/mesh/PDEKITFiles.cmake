list(APPEND PDEKIT_Mesh_HEADERS
Tria.hpp
CellGeometry.hpp
CellTransform.hpp
CellBuffer.hpp
DiscreteElemKey.hpp
DofCoordinates.hpp
CellMarker.hpp
ElementTopology.hpp
EntityDofRealign.hpp
EntityRealignCode.hpp
EntityStatus.hpp
KeyCache.hpp
MeshPredicates.hpp
MeshConstants.hpp
MeshEntity.hpp
MeshEntityIterator.hpp
MeshIndex.hpp
MeshConfig.hpp
Point.hpp
TopologyPredicates.hpp
)

list(APPEND PDEKIT_Mesh_SOURCES
CellTransform.cpp
CellMarker.cpp
DiscreteElemKey.cpp
EntityDofRealign.cpp
EntityRealignCode.cpp
EntityStatus.cpp
MeshPredicates.cpp
MeshEntity.cpp
MeshEntityIterator.cpp
TopologyPredicates.cpp
)
