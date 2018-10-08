#ifndef PDEKIT_Mesh_Internal_Trace_Tuple_Factory_hpp
#define PDEKIT_Mesh_Internal_Trace_Tuple_Factory_hpp

#include "common/FactoryPool.hpp"
#include "common/FactoryT.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"
#include "mesh/local_topology/TraceTupleBase.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

class TraceTupleFactorySetup
{
  public:
  static void setup_instance(
      common::FactoryT<TraceTupleBase, std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>>
          &instance);
};

} // namespace internal

typedef common::Singleton<
    common::FactoryT<internal::TraceTupleBase,
                     std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>>,
    internal::TraceTupleFactorySetup>
    TraceTupleFactory;

std::ostream &operator<<(
    std::ostream &os, const std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint> &tag_tuple);

} // namespace mesh

} // namespace pdekit

#endif
