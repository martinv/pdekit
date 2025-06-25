#include <iostream>

#include "mesh/EntityStatus.hpp"

namespace pdekit
{

namespace mesh
{

const std::string EntityStatusInfo::Names[5] = {"NotActive", "Active", "PendingRefinement",
                                                "PendingCoarsening", "PendingRemoval"};

std::ostream &operator<<(std::ostream &os, EntityStatus const status)
{
  os << EntityStatusInfo::Names[static_cast<std::underlying_type<EntityStatus>::type>(status)];
  return os;
}

} // namespace mesh

} // namespace pdekit
