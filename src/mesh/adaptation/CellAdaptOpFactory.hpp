#ifndef PDEKIT_Mesh_Adaptation_Cell_Adapt_Op_Factory_hpp
#define PDEKIT_Mesh_Adaptation_Cell_Adapt_Op_Factory_hpp

#include "common/FactoryPool.hpp"
#include "common/FactoryT.hpp"
#include "mesh/adaptation/CellAdaptOpBase.hpp"
#include "mesh/adaptation/CellAdaptOpTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

namespace internal
{

class CellAdaptOpFactorySetup
{
  public:
  static void setup_instance(common::FactoryT<CellAdaptOpBase, CellAdaptOpTag> &instance);
};

} // namespace internal

} // namespace adapt

typedef common::Singleton<common::FactoryT<adapt::internal::CellAdaptOpBase, adapt::CellAdaptOpTag>,
                          adapt::internal::CellAdaptOpFactorySetup>
    CellAdaptOpFactory;

} // namespace mesh

} // namespace pdekit

#endif
