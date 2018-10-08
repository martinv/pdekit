list(APPEND PDEKIT_Mesh_Containers_HEADERS
  containers/CellPath.hpp
  containers/TriaCells.hpp
  containers/TriaFacets.hpp
  containers/DofMap.hpp
  containers/MeshBoundary.hpp
  containers/MeshStatistics.hpp
  containers/MeshTopologyAdaptAlgorithm.hpp
  containers/TopologyAlgorithms.hpp
)

list(APPEND PDEKIT_Mesh_Containers_SOURCES
  containers/CellPath.cpp
  containers/MeshStatistics.cpp
  containers/TopologyAlgorithms.cpp
)
