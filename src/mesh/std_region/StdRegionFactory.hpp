#ifndef PDEKIT_Mesh_Std_Region_Factory_hpp
#define PDEKIT_Mesh_Std_Region_Factory_hpp

#include "common/FactoryT.hpp"
#include "common/Singleton.hpp"
#include "mesh/std_region/StdRegionBuilder.hpp"

namespace pdekit
{

namespace mesh
{

class StdRegionFactorySetup
{
  public:
  static void setup_instance(common::FactoryT<StdRegionBuilder, PointSetTag> &factory);
};

} // namespace mesh

typedef common::Singleton<common::FactoryT<mesh::StdRegionBuilder, mesh::PointSetTag>,
                          mesh::StdRegionFactorySetup>
    StdRegionFactory;

} // namespace pdekit

#endif
