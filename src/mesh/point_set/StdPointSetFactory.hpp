#ifndef PDEKIT_Mesh_Std_Point_Set_Factory_hpp
#define PDEKIT_Mesh_Std_Point_Set_Factory_hpp

#include "common/FactoryT.hpp"
#include "common/Singleton.hpp"
#include "mesh/point_set/StdPointSetBase.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

class StdPointSetFactorySetup
{
  public:
  static void setup_instance(common::FactoryT<StdPointSetBase, mesh::PointSetTag> &instance);
};

} // namespace mesh

typedef common::Singleton<common::FactoryT<mesh::StdPointSetBase, mesh::PointSetTag>,
                          mesh::StdPointSetFactorySetup>
    StdPointSetFactory;

} // namespace pdekit

#endif
