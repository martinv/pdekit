list(APPEND PDEKIT_Mesh_Iterators_HEADERS
  iterators/CellDofIterator.hpp
  iterators/CellTopologyIteratorFilter.hpp
  iterators/DofIteratorFilter.hpp
  iterators/BdryDofIterator.hpp
  iterators/CellTopologyIterator.hpp
  iterators/IteratorInterface.hpp
  iterators/TraceTopologyIterator.hpp
)

list(APPEND PDEKIT_Mesh_Iterators_SOURCES
  iterators/CellTopologyIteratorFilter.cpp
  iterators/DofIteratorFilter.cpp
)
