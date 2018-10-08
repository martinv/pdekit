#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"
#include "common/FactoryPool.hpp"
#include "mesh/point_set/QuadratureLineAdaptTransformAlgo.hpp"

namespace pdekit
{

namespace mesh
{

void QuadAdaptTransformAlgoFactSetup::setup_instance(
    common::FactoryT<QuadratureAdaptTransformAlgoBase, mesh::CellTransform> &instance)
{

  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> qd_adapt_trans_factory =
      QuadratureAdaptTransformAlgoFactory::shared_ptr_to_instance();
  fpool.register_factory(QuadratureAdaptTransformAlgoBase::type_name(), qd_adapt_trans_factory);

  // Fill the factory
  instance.register_builder<QuadLineToTwoSegmentsTransformAlgo>(
      mesh::CellTransform::UNIFORM_REFINE);
}

} // namespace mesh

template class common::Singleton<
    common::FactoryT<mesh::QuadratureAdaptTransformAlgoBase, mesh::CellTransform>,
    mesh::QuadAdaptTransformAlgoFactSetup>;

} // namespace pdekit
