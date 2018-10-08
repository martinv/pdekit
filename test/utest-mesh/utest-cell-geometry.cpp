/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_geometry_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "mesh/CellGeometry.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_geometry_utest)
{
  /// Node ids of a P2 triangle
  std::vector<Uint> p2_tri_node_ids_raw = {0, 1, 2, 3, 4, 5};
  common::ArrayView<const Uint, _1D, Uint> p2_tri_node_ids_view(p2_tri_node_ids_raw.data(),
                                                                p2_tri_node_ids_raw.size());
  const mesh::PointSetTag p2_tri_tag(ElemShape::Triag, P2, PointSetID::Equidist);

  const mesh::MeshEntity p2_tri(p2_tri_node_ids_view, 0, p2_tri_tag);

  std::vector<Real> p2_tri_node_coord_raw = {-1., -1., 1., -1., -1., 1., 0., -1., 0., 0., -1., 0.};
  const math::DenseVecView<const Real> p2_tri_node_coord_view(p2_tri_node_coord_raw.data(),
                                                              p2_tri_node_coord_raw.size());

  const mesh::MeshEntity p2_tri_edge = p2_tri.sub_entity(_1D, 2);

  const mesh::CellGeometry<_2D> p2_tri_cell_geo(p2_tri_node_coord_view, p2_tri_edge);

  const mesh::CellGeometry<_2D>::node_view_t node = p2_tri_cell_geo.const_node_view(2);
  std::cout << "One node: " << node << std::endl;

  std::cout << p2_tri_cell_geo << std::endl;
}

// ----------------------------------------------------------------------------
