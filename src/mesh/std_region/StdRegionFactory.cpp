#include "mesh/std_region/StdRegionFactory.hpp"
#include "common/FactoryPool.hpp"
#include "common/Provider.hpp"
#include "mesh/std_region/EquidistStdRegion.hpp"
#include "mesh/std_region/StdRegionBuilderT.hpp"
#include "mesh/std_region/UndefinedStdRegion.hpp"
#include "mesh/std_region/WarpblendStdRegion.hpp"

namespace pdekit
{

// template class
// common::AddToFactory<common::FactoryT<mesh::sf::ShapeFunction>,mesh::sf::LagrangeP2Triag2D>;

namespace mesh
{

/*
common::Provider<StdRegionFactory, StdRegionT<UndefinedStdRegion>>
undef_entity_provider(PointSetTag(ElemShapeID::Undefined, P0, Undefined));

common::Provider<StdRegionFactory, StdRegionT<EquidistStdRegion<Line, P1>>>
p1_equidist_line_provider(PointSetTag(Line, P1, Equidist));
*/

void StdRegionFactorySetup::setup_instance(common::FactoryT<StdRegionBuilder, PointSetTag> &factory)
{
  //  std::cout << "SF setup" << std::endl;

  // Register itself in FactoryPool:
  common::FactoryPool::instance_type &fpool = common::FactoryPool::instance();
  std::shared_ptr<common::AbstractFactoryBase> std_region_factory =
      StdRegionFactory::shared_ptr_to_instance();
  fpool.register_factory(StdRegionBuilder::type_name(), std_region_factory);

  // Fill the factory
  factory.register_builder<StdRegionBuilderT<UndefinedStdRegion>, const Uint>(
      PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), P0);

  // --------------------------------------------------------------------------
  // LINES
  // --------------------------------------------------------------------------

  // Standard region for lines, equidistant node distribution
  for (Uint p = P1; p <= P15; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionLine>, const Uint>(
        PointSetTag(ElemShape::Line, p, PointSetID::Equidist), p);
  }

  // Standard region for lines, warpblend node distribution
  for (Uint p = P1; p <= P15; ++p)
  {
    factory.register_builder<StdRegionBuilderT<WarpblendStdRegionLine>, const Uint>(
        PointSetTag(ElemShape::Line, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // TRIANGLES
  // --------------------------------------------------------------------------

  // Standard region for triangles, equidistant node distribution
  for (Uint p = P1; p <= P15; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionTriag>, const Uint>(
        PointSetTag(ElemShape::Triag, p, PointSetID::Equidist), p);
  }

  // Standard region for triangles, warpblend node distribution
  for (Uint p = P1; p <= P15; ++p)
  {
    factory.register_builder<StdRegionBuilderT<WarpblendStdRegionTriag>, const Uint>(
        PointSetTag(ElemShape::Triag, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // QUADRILATERALS
  // --------------------------------------------------------------------------

  // Standard region for quadrilaterals, equidistant node distribution
  for (Uint p = P1; p <= P10; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionQuad>, const Uint>(
        PointSetTag(ElemShape::Quad, p, PointSetID::Equidist), p);
  }

  // Standard region for quadrilaterals, warpblend node distribution
  for (Uint p = P1; p <= P10; ++p)
  {
    factory.register_builder<StdRegionBuilderT<WarpblendStdRegionQuad>, const Uint>(
        PointSetTag(ElemShape::Quad, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // TETRAHEDRA
  // --------------------------------------------------------------------------

  // Standard region for tetrahedra, equidistant node distribution
  for (Uint p = P1; p <= P9; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionTetra>, const Uint>(
        PointSetTag(ElemShape::Tetra, p, PointSetID::Equidist), p);
  }

  // Standard region for tetrahedra, warpblend node distribution
  for (Uint p = P1; p <= P9; ++p)
  {
    factory.register_builder<StdRegionBuilderT<WarpblendStdRegionTetra>, const Uint>(
        PointSetTag(ElemShape::Tetra, p, PointSetID::Warpblend), p);
  }

  // --------------------------------------------------------------------------
  // HEXAHEDRA
  // --------------------------------------------------------------------------

  // Standard region for hexahedra, equidistant node distribution
  for (Uint p = P1; p <= P6; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionHexa>, const Uint>(
        PointSetTag(ElemShape::Hexa, p, PointSetID::Equidist), p);
  }

  // --------------------------------------------------------------------------
  // PYRAMIDS
  // --------------------------------------------------------------------------

  // Standard region for pyramids, equidistant node distribution
  for (Uint p = P1; p <= P6; ++p)
  {
    factory.register_builder<StdRegionBuilderT<EquidistStdRegionPyramid>, const Uint>(
        PointSetTag(ElemShape::Pyramid, p, PointSetID::Equidist), p);
  }
}

} // namespace mesh

// Force instantantiation of some classes
template class common::FactoryT<mesh::StdRegionBuilder, mesh::PointSetTag>;
template class common::Singleton<common::FactoryT<mesh::StdRegionBuilder, mesh::PointSetTag>,
                                 mesh::StdRegionFactorySetup>;

} // namespace pdekit
