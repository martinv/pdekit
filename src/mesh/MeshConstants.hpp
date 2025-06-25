#ifndef PDEKIT_Mesh_Mesh_Constants_hpp
#define PDEKIT_Mesh_Mesh_Constants_hpp

#include "common/PDEKit.hpp"
#include <limits>

namespace pdekit
{

namespace mesh
{

enum DofStorageType
{
  UndefinedDofStorage         = 0,
  CellCentered                = 1,
  VertexCenteredContinuous    = 2,
  VertexCenteredDiscontinuous = 3
};

/*
const Uint INVALID_CELL_ID = std::numeric_limits<Uint>::max();
const Uint INVALID_NODE_ID = std::numeric_limits<Uint>::max() - 1;
const Uint INVALID_ENTITY_ID = std::numeric_limits<Uint>::max() - 2;
const Uint INVALID_DOF_ID = std::numeric_limits<Uint>::max() - 3;
const Uint INVALID_REF_ENTITY_ID = std::numeric_limits<SUint>::max();
*/

// Assigning std::numeric_limits<>::max() to enum works in C++11 because
// std::numeric_limits<>::max() return a constexpr T here

enum
{
  INVALID_CELL_ID       = std::numeric_limits<Uint>::max(),
  INVALID_FACET_ID      = std::numeric_limits<Uint>::max() - 1,
  INVALID_EDGE_ID       = std::numeric_limits<Uint>::max() - 2,
  INVALID_NODE_ID       = std::numeric_limits<Uint>::max() - 3,
  INVALID_ENTITY_ID     = std::numeric_limits<Uint>::max() - 4,
  INVALID_DOF_ID        = std::numeric_limits<Uint>::max() - 5,
  INVALID_REF_ENTITY_ID = std::numeric_limits<SUint>::max(),
  INVALID_LOC_ENTITY_ID = std::numeric_limits<SUint>::max() - 1
};

} // namespace mesh

} // namespace pdekit

#endif
