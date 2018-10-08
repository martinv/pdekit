#include "mesh/point_set/StdPointSetFactory.hpp"
#include "common/FactoryPool.hpp"
#include "mesh/point_set/DofPointSetHexaEquidist.hpp"
#include "mesh/point_set/DofPointSetLineEquidist.hpp"
#include "mesh/point_set/DofPointSetLineWarpblend.hpp"
#include "mesh/point_set/DofPointSetPyramidEquidist.hpp"
#include "mesh/point_set/DofPointSetQuadEquidist.hpp"
#include "mesh/point_set/DofPointSetQuadWarpblend.hpp"
#include "mesh/point_set/DofPointSetTetraEquidist.hpp"
#include "mesh/point_set/DofPointSetTetraWarpblend.hpp"
#include "mesh/point_set/DofPointSetTriagEquidist.hpp"
#include "mesh/point_set/DofPointSetTriagWarpblend.hpp"
#include "mesh/point_set/QdPointSetHexaGauss.hpp"
#include "mesh/point_set/QdPointSetLineGauss.hpp"
#include "mesh/point_set/QdPointSetLineGaussLobatto.hpp"
#include "mesh/point_set/QdPointSetQuadFaceGauss.hpp"
#include "mesh/point_set/QdPointSetQuadGauss.hpp"
#include "mesh/point_set/QdPointSetTetraFaceGauss.hpp"
#include "mesh/point_set/QdPointSetTetraGauss.hpp"
#include "mesh/point_set/QdPointSetTriagFaceGauss.hpp"
#include "mesh/point_set/QdPointSetTriagGauss.hpp"

namespace pdekit
{

namespace mesh
{

void StdPointSetFactorySetup::setup_instance(
    common::FactoryT<StdPointSetBase, mesh::PointSetTag> &instance)
{
  //  std::cout << "Quadrature factory setup" << std::endl;

  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> qd_factory =
      StdPointSetFactory::shared_ptr_to_instance();
  fpool.register_factory(StdPointSetBase::type_name(), qd_factory);

  // Fill the factory

  // --------------------------------------------------------------------------
  // LINES
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P15; ++p)
  {
    instance.register_builder<DofPointSetLineEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Line, p, PointSetID::Equidist), p);
  }

  for (Uint p = P1; p <= P15; ++p)
  {
    instance.register_builder<DofPointSetLineWarpblend, const Uint>(
        mesh::PointSetTag(ElemShape::Line, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // TRIANGLES
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P15; ++p)
  {
    instance.register_builder<DofPointSetTriagEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Triag, p, PointSetID::Equidist), p);
  }

  for (Uint p = P1; p <= P15; ++p)
  {
    instance.register_builder<DofPointSetTriagWarpblend, const Uint>(
        mesh::PointSetTag(ElemShape::Triag, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // QUADRILATERALS
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P10; ++p)
  {
    instance.register_builder<DofPointSetQuadEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Quad, p, PointSetID::Equidist), p);
  }

  for (Uint p = P1; p <= P10; ++p)
  {
    instance.register_builder<DofPointSetQuadWarpblend, const Uint>(
        mesh::PointSetTag(ElemShape::Quad, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // TETRAHEDRA
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P9; ++p)
  {
    instance.register_builder<DofPointSetTetraEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Tetra, p, PointSetID::Equidist), p);
  }

  for (Uint p = P1; p <= P9; ++p)
  {
    instance.register_builder<DofPointSetTetraWarpblend, const Uint>(
        mesh::PointSetTag(ElemShape::Tetra, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // HEXAHEDRA
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P6; ++p)
  {
    instance.register_builder<DofPointSetHexaEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Hexa, p, PointSetID::Equidist), p);
  }

  // --------------------------------------------------------------------------
  // PYRAMIDS
  // --------------------------------------------------------------------------

  for (Uint p = P1; p <= P6; ++p)
  {
    instance.register_builder<DofPointSetPyramidEquidist, const Uint>(
        mesh::PointSetTag(ElemShape::Pyramid, p, PointSetID::Equidist), p);
  }

  instance.register_builder<QdPointSetP1LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P1, PointSetID::Gauss));
  instance.register_builder<QdPointSetP2LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P2, PointSetID::Gauss));
  instance.register_builder<QdPointSetP3LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P3, PointSetID::Gauss));
  instance.register_builder<QdPointSetP4LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P4, PointSetID::Gauss));
  instance.register_builder<QdPointSetP5LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P5, PointSetID::Gauss));
  instance.register_builder<QdPointSetP6LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P6, PointSetID::Gauss));
  instance.register_builder<QdPointSetP7LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P7, PointSetID::Gauss));
  instance.register_builder<QdPointSetP8LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P8, PointSetID::Gauss));
  instance.register_builder<QdPointSetP9LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P9, PointSetID::Gauss));
  instance.register_builder<QdPointSetP10LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P10, PointSetID::Gauss));
  instance.register_builder<QdPointSetP11LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P11, PointSetID::Gauss));
  instance.register_builder<QdPointSetP12LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P12, PointSetID::Gauss));
  instance.register_builder<QdPointSetP13LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P13, PointSetID::Gauss));
  instance.register_builder<QdPointSetP14LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P14, PointSetID::Gauss));
  instance.register_builder<QdPointSetP15LineGauss>(
      mesh::PointSetTag(ElemShape::Line, P15, PointSetID::Gauss));

  instance.register_builder<QdPointSetP1LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P1, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP2LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P2, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP3LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P3, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP4LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P4, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP5LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P5, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP6LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P6, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP7LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P7, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP8LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P8, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP9LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P9, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP10LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P10, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP11LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P11, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP12LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P12, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP13LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P13, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP14LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P14, PointSetID::GaussLobatto));
  instance.register_builder<QdPointSetP15LineGaussLobatto>(
      mesh::PointSetTag(ElemShape::Line, P15, PointSetID::GaussLobatto));

  instance.register_builder<QdPointSetP1TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Gauss));
  instance.register_builder<QdPointSetP2TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P2, PointSetID::Gauss));
  instance.register_builder<QdPointSetP3TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P3, PointSetID::Gauss));
  instance.register_builder<QdPointSetP4TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P4, PointSetID::Gauss));
  instance.register_builder<QdPointSetP5TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P5, PointSetID::Gauss));
  instance.register_builder<QdPointSetP6TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P6, PointSetID::Gauss));
  instance.register_builder<QdPointSetP7TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P7, PointSetID::Gauss));
  instance.register_builder<QdPointSetP8TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P8, PointSetID::Gauss));
  instance.register_builder<QdPointSetP9TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P9, PointSetID::Gauss));
  instance.register_builder<QdPointSetP10TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P10, PointSetID::Gauss));
  instance.register_builder<QdPointSetP11TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P11, PointSetID::Gauss));
  instance.register_builder<QdPointSetP12TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P12, PointSetID::Gauss));
  instance.register_builder<QdPointSetP13TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P13, PointSetID::Gauss));
  instance.register_builder<QdPointSetP14TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P14, PointSetID::Gauss));
  instance.register_builder<QdPointSetP15TriagGauss>(
      mesh::PointSetTag(ElemShape::Triag, P15, PointSetID::Gauss));

  // Quadrature on triangle faces depends on line quadrature,
  // so only those quadrature rules whose order is NOT HIGHER
  // THAN THE ORDER OF LINE QUADRATURE RULES can be registered
  instance.register_builder<QdPointSetTriagFaceGauss<P1>>(
      mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P2>>(
      mesh::PointSetTag(ElemShape::Triag, P2, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P3>>(
      mesh::PointSetTag(ElemShape::Triag, P3, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P4>>(
      mesh::PointSetTag(ElemShape::Triag, P4, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P5>>(
      mesh::PointSetTag(ElemShape::Triag, P5, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P6>>(
      mesh::PointSetTag(ElemShape::Triag, P6, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P7>>(
      mesh::PointSetTag(ElemShape::Triag, P7, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P8>>(
      mesh::PointSetTag(ElemShape::Triag, P8, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P9>>(
      mesh::PointSetTag(ElemShape::Triag, P9, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P10>>(
      mesh::PointSetTag(ElemShape::Triag, P10, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P11>>(
      mesh::PointSetTag(ElemShape::Triag, P11, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P12>>(
      mesh::PointSetTag(ElemShape::Triag, P12, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P13>>(
      mesh::PointSetTag(ElemShape::Triag, P13, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P14>>(
      mesh::PointSetTag(ElemShape::Triag, P14, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTriagFaceGauss<P15>>(
      mesh::PointSetTag(ElemShape::Triag, P15, PointSetID::FaceGauss));

  instance.register_builder<QdPointSetP1QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P1, PointSetID::Gauss));
  instance.register_builder<QdPointSetP2QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P2, PointSetID::Gauss));
  instance.register_builder<QdPointSetP3QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P3, PointSetID::Gauss));
  instance.register_builder<QdPointSetP4QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P4, PointSetID::Gauss));
  instance.register_builder<QdPointSetP5QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P5, PointSetID::Gauss));
  instance.register_builder<QdPointSetP6QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P6, PointSetID::Gauss));
  instance.register_builder<QdPointSetP7QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P7, PointSetID::Gauss));
  instance.register_builder<QdPointSetP8QuadGauss>(
      mesh::PointSetTag(ElemShape::Quad, P8, PointSetID::Gauss));

  // Quadrature on faces of reference quadrilateral depends on line
  // quadrature, so only those quadrature rules whose order is NOT
  // HIGHER THAN THE ORDER OF LINE QUADRATURE RULES can be registered
  instance.register_builder<QdPointSetQuadFaceGauss<P1>>(
      mesh::PointSetTag(ElemShape::Quad, P1, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P2>>(
      mesh::PointSetTag(ElemShape::Quad, P2, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P3>>(
      mesh::PointSetTag(ElemShape::Quad, P3, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P4>>(
      mesh::PointSetTag(ElemShape::Quad, P4, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P5>>(
      mesh::PointSetTag(ElemShape::Quad, P5, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P6>>(
      mesh::PointSetTag(ElemShape::Quad, P6, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P7>>(
      mesh::PointSetTag(ElemShape::Quad, P7, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P8>>(
      mesh::PointSetTag(ElemShape::Quad, P8, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P9>>(
      mesh::PointSetTag(ElemShape::Quad, P9, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetQuadFaceGauss<P10>>(
      mesh::PointSetTag(ElemShape::Quad, P10, PointSetID::FaceGauss));

  instance.register_builder<QdPointSetP1TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P1, PointSetID::Gauss));
  instance.register_builder<QdPointSetP2TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P2, PointSetID::Gauss));
  instance.register_builder<QdPointSetP3TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P3, PointSetID::Gauss));
  instance.register_builder<QdPointSetP4TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P4, PointSetID::Gauss));
  instance.register_builder<QdPointSetP5TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P5, PointSetID::Gauss));
  instance.register_builder<QdPointSetP6TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P6, PointSetID::Gauss));
  instance.register_builder<QdPointSetP7TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P7, PointSetID::Gauss));
  instance.register_builder<QdPointSetP8TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P8, PointSetID::Gauss));
  instance.register_builder<QdPointSetP9TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P9, PointSetID::Gauss));
  instance.register_builder<QdPointSetP10TetraGauss>(
      mesh::PointSetTag(ElemShape::Tetra, P10, PointSetID::Gauss));

  instance.register_builder<QdPointSetTetraFaceGauss<P1>>(
      mesh::PointSetTag(ElemShape::Tetra, P1, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P2>>(
      mesh::PointSetTag(ElemShape::Tetra, P2, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P3>>(
      mesh::PointSetTag(ElemShape::Tetra, P3, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P4>>(
      mesh::PointSetTag(ElemShape::Tetra, P4, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P5>>(
      mesh::PointSetTag(ElemShape::Tetra, P5, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P6>>(
      mesh::PointSetTag(ElemShape::Tetra, P6, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P7>>(
      mesh::PointSetTag(ElemShape::Tetra, P7, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P8>>(
      mesh::PointSetTag(ElemShape::Tetra, P8, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P9>>(
      mesh::PointSetTag(ElemShape::Tetra, P9, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P10>>(
      mesh::PointSetTag(ElemShape::Tetra, P10, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P11>>(
      mesh::PointSetTag(ElemShape::Tetra, P11, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P12>>(
      mesh::PointSetTag(ElemShape::Tetra, P12, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P13>>(
      mesh::PointSetTag(ElemShape::Tetra, P13, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P14>>(
      mesh::PointSetTag(ElemShape::Tetra, P14, PointSetID::FaceGauss));
  instance.register_builder<QdPointSetTetraFaceGauss<P15>>(
      mesh::PointSetTag(ElemShape::Tetra, P15, PointSetID::FaceGauss));

  instance.register_builder<QdPointSetP1HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P1, PointSetID::Gauss));
  instance.register_builder<QdPointSetP2HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P2, PointSetID::Gauss));
  instance.register_builder<QdPointSetP3HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P3, PointSetID::Gauss));
  instance.register_builder<QdPointSetP4HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P4, PointSetID::Gauss));
  instance.register_builder<QdPointSetP5HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P5, PointSetID::Gauss));
  instance.register_builder<QdPointSetP6HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P6, PointSetID::Gauss));
  instance.register_builder<QdPointSetP7HexaGauss>(
      mesh::PointSetTag(ElemShape::Hexa, P7, PointSetID::Gauss));
}

} // namespace mesh

template class common::FactoryT<mesh::StdPointSetBase>;
template class common::Singleton<common::FactoryT<mesh::StdPointSetBase, mesh::PointSetTag>,
                                 mesh::StdPointSetFactorySetup>;

} // namespace pdekit
