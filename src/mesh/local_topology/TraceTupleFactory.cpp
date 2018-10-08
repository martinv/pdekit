#include "mesh/local_topology/TraceTupleFactory.hpp"
#include "common/FactoryPool.hpp"
#include "mesh/local_topology/TraceTupleLine.hpp"
#include "mesh/local_topology/TraceTupleQuad.hpp"
#include "mesh/local_topology/TraceTupleTriag.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

void TraceTupleFactorySetup::setup_instance(
    common::FactoryT<TraceTupleBase, std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>>
        &instance)
{
  //  std::cout << "SF setup" << std::endl;

  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> facet_incidence_factory =
      TraceTupleFactory::shared_ptr_to_instance();
  fpool.register_factory(TraceTupleBase::type_name(), facet_incidence_factory);

  // IMPORTANT:
  // Note that in the factory, the two CellSubdomainTags define what is the
  // mutual position of entities BEFORE adaptation

  instance.register_builder<LineL_L0_LineR_L0>(
      std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>(
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 0,
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 0));

  instance.register_builder<LineL_L1_LineR_L0PartType1Pos0>(
      std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>(
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 1,
          CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line), 0));

  instance.register_builder<LineL_L1_LineR_L0PartType1Pos1>(
      std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>(
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 1,
          CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line), 0));

  instance.register_builder<LineL_L0PartType1Pos0_LineR_L1>(
      std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>(
          CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 0, ElemShape::Line), 0,
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 1));

  instance.register_builder<LineL_L0PartType1Pos1_LineR_L1>(
      std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint>(
          CellSubdomainTag(ElemShape::Line, CellTransform::UNIFORM_REFINE, 1, ElemShape::Line), 0,
          CellSubdomainTag(ElemShape::Line, CellTransform::NO_TRANS, 0, ElemShape::Line), 1));
}

} // namespace internal

std::ostream &operator<<(
    std::ostream &os, const std::tuple<CellSubdomainTag, Uint, CellSubdomainTag, Uint> &tag_tuple)
{
  os << std::get<0>(tag_tuple).as_string() << "_L" << std::get<1>(tag_tuple) << " <-> "
     << std::get<2>(tag_tuple).as_string() << "_L" << std::get<3>(tag_tuple);
  return os;
}

} // namespace mesh

// Force instantantiation of certain classes
// This has to be done in a namespace that is enclosing common:: and mesh:: (?)
template class common::FactoryT<
    mesh::internal::TraceTupleBase,
    std::tuple<mesh::CellSubdomainTag, Uint, mesh::CellSubdomainTag, Uint>>;
template class common::Singleton<
    common::FactoryT<mesh::internal::TraceTupleBase,
                     std::tuple<mesh::CellSubdomainTag, Uint, mesh::CellSubdomainTag, Uint>>,
    mesh::internal::TraceTupleFactorySetup>;

} // namespace pdekit
