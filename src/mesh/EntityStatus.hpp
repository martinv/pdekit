#ifndef PDEKIT_Mesh_Entity_Status_hpp
#define PDEKIT_Mesh_Entity_Status_hpp

#include <iosfwd>

namespace pdekit
{

namespace mesh
{

enum class EntityStatus : unsigned short
{
  NotActive         = 0,
  Active            = 1,
  PendingRefinement = 2,
  PendingCoarsening = 3,
  PendingRemoval    = 4
};

struct EntityStatusInfo
{
  static const std::string Names[5];
};

std::ostream &operator<<(std::ostream &os, EntityStatus const status);

} // namespace mesh

} // namespace pdekit

#endif
