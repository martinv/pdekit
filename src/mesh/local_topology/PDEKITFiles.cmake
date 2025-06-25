list(APPEND PDEKIT_Mesh_Local_Topology_HEADERS
  local_topology/TraceTupleBase.hpp
  local_topology/TraceTupleLine.hpp
  local_topology/TraceTupleTriag.hpp
  local_topology/TraceTupleQuad.hpp
  local_topology/TraceTupleFactory.hpp
  local_topology/TraceEntityTuple.hpp
  local_topology/CellSubdomainTag.hpp
  local_topology/TraceIncidences.hpp
)

list(APPEND PDEKIT_Mesh_Local_Topology_SOURCES
  local_topology/TraceTupleBase.cpp
  local_topology/TraceTupleLine.cpp
  local_topology/TraceTupleTriag.cpp
  local_topology/TraceTupleQuad.cpp
  local_topology/TraceTupleFactory.cpp
  local_topology/TraceEntityTuple.cpp
  local_topology/CellSubdomainTag.cpp
  local_topology/TraceIncidences.cpp
)
