#ifndef PDEKIT_Mesh_Quadrature_Adapt_Transform_Algo_Factory_hpp
#define PDEKIT_Mesh_Quadrature_Adapt_Transform_Algo_Factory_hpp

#include "common/FactoryT.hpp"
#include "common/Singleton.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoBase.hpp"

namespace pdekit
{

namespace mesh
{

class QuadAdaptTransformAlgoFactSetup
{
  public:
  static void setup_instance(
      common::FactoryT<QuadratureAdaptTransformAlgoBase, mesh::CellTransform> &instance);
};

} // namespace mesh

typedef common::Singleton<
    common::FactoryT<mesh::QuadratureAdaptTransformAlgoBase, mesh::CellTransform>,
    mesh::QuadAdaptTransformAlgoFactSetup>
    QuadratureAdaptTransformAlgoFactory;

} // namespace pdekit

#endif
