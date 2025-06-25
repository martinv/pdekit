/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE polynomial_metric_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <iomanip>
#include <iostream>
#include <memory>

#include "interpolation/PolyCache.hpp"
#include "interpolation/PolyMetric.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------
// This function is supposed to be used with quadratures that place at
// least some of their points in element interior, such as standard
// Gauss quadrature
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void test_metric(typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
                 const PointSetID quadrature_type, const Real reference_volume_value,
                 const Real tol)
{
  std::cout << std::endl;
  std::cout << "**********************************************************" << std::endl;
  std::cout << " Test metric [" << MeshConfig::TDIM << "D, " << DIM << "D]" << std::endl;
  std::cout << "**********************************************************" << std::endl;

  clock_t t1, t2, t3;
  Real elapsed;

  const auto geo_sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto sol_sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order + 1, ModalBasis::Modal);
  };

  const auto eval_pt_set_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, P4, quadrature_type);
  };

  interpolation::FunctionSpace<MeshConfig, DIM> fs_geometry;
  fs_geometry.set_reference_fe_values(common::make_iter_range(cell_dofs.cbegin(), cell_dofs.cend()),
                                      geo_sf_generator, eval_pt_set_generator);

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> geometry_metric;

  geometry_cache.allocate(fs_geometry.discrete_elements().cbegin(),
                          fs_geometry.discrete_elements().cend(), cell_dofs.nb_active_cells());
  // For the moment, we use the same interpolation space for geometry and
  // solution
  geometry_metric.allocate_buffer(fs_geometry.discrete_elements().cbegin(),
                                  fs_geometry.discrete_elements().cend(),
                                  cell_dofs.nb_active_cells());

  BOOST_CHECK_EQUAL(geometry_metric.max_nb_blocks_in_buffer(), cell_dofs.nb_active_cells());

  interpolation::PolyCacheDescription pdesc;
  interpolation::PolyCache poly_cache;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell         = cell_dofs.active_cell(mesh::ActiveIdx(c));
    const mesh::PointSetTag std_reg_tag = cell.pt_set_id();
    const ElemShape cell_shape          = std_reg_tag.elem_shape();
    const Uint cell_order               = std_reg_tag.poly_order();
    const mesh::PointSetTagExt std_reg_tag_ext(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::sf::SFTag sol_basis_tag = sol_sf_generator(cell_shape, cell_order);

    const mesh::PointSetTag eval_pt_set_tag = eval_pt_set_generator(cell_shape, cell_order);
    const mesh::PointSetTagExt eval_pt_set_tag_ext(eval_pt_set_tag, P0,
                                                   mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey poly_key(std_reg_tag_ext, sol_basis_tag, eval_pt_set_tag_ext);
    pdesc.add_to_block_count(poly_key, 1);
  }

  poly_cache.allocate(pdesc);

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const mesh::CellGeometry<MeshConfig::GDIM> cell_coords = tcell_view.coordinates();

    const mesh::PointSetTag std_reg_tag = cell.pt_set_id();
    const ElemShape cell_shape          = std_reg_tag.elem_shape();
    const Uint cell_order               = std_reg_tag.poly_order();
    const mesh::PointSetTagExt std_reg_tag_ext =
        mesh::PointSetTagExt(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::sf::SFTag geo_basis_tag = geo_sf_generator(cell_shape, cell_order);
    const mesh::sf::SFTag sol_basis_tag = sol_sf_generator(cell_shape, cell_order);

    const mesh::PointSetTag eval_pt_set_tag = eval_pt_set_generator(cell_shape, cell_order);
    const mesh::PointSetTagExt eval_pt_set_tag_ext =
        mesh::PointSetTagExt(eval_pt_set_tag, P0, mesh::CellTransform::NO_TRANS, 0u);

    const mesh::DiscreteElemKey geo_key(std_reg_tag_ext, geo_basis_tag, eval_pt_set_tag_ext);
    const mesh::DiscreteElemKey poly_key(std_reg_tag_ext, sol_basis_tag, eval_pt_set_tag_ext);

    geometry_cache.push_back_to_buffer(cell_coords, geo_key);

    poly_cache.push_back_to_buffer(poly_key);
  }

  t1 = clock();

  geometry_metric.evaluate(geometry_cache, interpolation::RebuildMetricIndex{true});

  t2 = clock();

  using geometry_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             DIM>::cellwise_metric;

  math::DenseSVec<Real, MeshConfig::GDIM> dsf_ref;
  math::DenseSVec<Real, MeshConfig::GDIM> dsf_phys;
  math::DenseDMat<Real> check_sums;

  std::vector<std::reference_wrapper<const math::DenseDMat<Real>>> dV;
  //  dV.resize(MeshConfig::GDIM);

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const interpolation::FEValues &fe_values = poly_cache.fe_values(c);
    // math::DenseDMat<Real> const &V = fe_values.Vandermonde();

    for (Uint d = 0; d < MeshConfig::GDIM; ++d)
    {
      dV.push_back(std::cref(fe_values.deriv_Vandermonde(d)));
    }

    /*
    for (Uint i = 0; i < dV.size(); ++i)
    {
      std::cout << dV[i].get() << std::endl;
    }
    */

    interpolation::PolyMetric<MeshConfig, DIM> poly_metric;

    check_sums.resize(fe_values.nb_nodes(), MeshConfig::GDIM);

    geometry_metric_type geo_met = geometry_metric.cellwise_values(c);
  }

  t3 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + computation of jacobians) = " << elapsed << " s"
            << std::endl;
  elapsed = ((Real)(t3 - t2)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + reuse of jacobians) = " << elapsed << " s" << std::endl;

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);
  // std::cout << "Domain measure = " << std::setw(10) << domain_volume <<
  // std::endl;

  // BOOST_CHECK_CLOSE(domain_volume, reference_volume_value, tol);

  std::cout << "**********************************************************" << std::endl;
}

// ----------------------------------------------------------------------------

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(polynomial_cache_description_utest)
{
  interpolation::PolyCacheDescription pdesc;

  BOOST_CHECK_EQUAL(pdesc.nb_blocks(), 0u);

  const mesh::PointSetTag elem_tag1(ElemShape::Quad, P2, PointSetID::Equidist);
  const mesh::PointSetTagExt elem_tag1_ext(elem_tag1, P0, mesh::CellTransform::NO_TRANS, 0u);
  const mesh::PointSetTag quad_tag1(ElemShape::Quad, P4, PointSetID::Gauss);
  const mesh::PointSetTagExt quad_tag1_ext(quad_tag1, P0, mesh::CellTransform::NO_TRANS, 0u);

  const mesh::PointSetTag elem_tag2(ElemShape::Triag, P3, PointSetID::Warpblend);
  const mesh::PointSetTagExt elem_tag2_ext(elem_tag2, P0, mesh::CellTransform::NO_TRANS, 0u);
  const mesh::PointSetTag quad_tag2(ElemShape::Triag, P6, PointSetID::Gauss);
  const mesh::PointSetTagExt quad_tag2_ext(quad_tag2, P0, mesh::CellTransform::NO_TRANS, 0u);

  const mesh::sf::SFTag basis_tag1(ElemShape::Quad, SFunc::Lagrange, P2, ModalBasis::Modal);
  const mesh::sf::SFTag basis_tag2a(ElemShape::Triag, SFunc::Lagrange, P3, ModalBasis::Modal);
  const mesh::sf::SFTag basis_tag2b(ElemShape::Triag, SFunc::Lagrange, P4, ModalBasis::Modal);

  for (Uint i = 0; i < 4; ++i)
  {
    pdesc.add_to_block_count(mesh::DiscreteElemKey(elem_tag1_ext, basis_tag1, quad_tag1_ext));
  }

  for (Uint i = 0; i < 4; ++i)
  {
    pdesc.add_to_block_count(mesh::DiscreteElemKey(elem_tag2_ext, basis_tag2a, quad_tag2_ext));
  }

  for (Uint i = 0; i < 2; ++i)
  {
    pdesc.add_to_block_count(mesh::DiscreteElemKey(elem_tag2_ext, basis_tag2b, quad_tag2_ext));
  }

  BOOST_CHECK_EQUAL(pdesc.nb_blocks(), 3u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(polynomial_cache_utest)
{
  const mesh::PointSetTag elem_tag1(ElemShape::Quad, P2, PointSetID::Equidist);
  const mesh::PointSetTagExt elem_tag1_ext(elem_tag1, P0, mesh::CellTransform::NO_TRANS, 0u);
  const mesh::PointSetTag quad_tag1(ElemShape::Quad, P4, PointSetID::Gauss);
  const mesh::PointSetTagExt quad_tag1_ext(quad_tag1, P0, mesh::CellTransform::NO_TRANS, 0u);

  const mesh::PointSetTag elem_tag2(ElemShape::Triag, P3, PointSetID::Warpblend);
  const mesh::PointSetTagExt elem_tag2_ext(elem_tag2, P0, mesh::CellTransform::NO_TRANS, 0u);
  const mesh::PointSetTag quad_tag2(ElemShape::Triag, P6, PointSetID::Gauss);
  const mesh::PointSetTagExt quad_tag2_ext(quad_tag2, P0, mesh::CellTransform::NO_TRANS, 0u);

  const mesh::sf::SFTag basis_tag1(ElemShape::Quad, SFunc::Lagrange, P2, ModalBasis::Modal);
  const mesh::sf::SFTag basis_tag2a(ElemShape::Triag, SFunc::Lagrange, P3, ModalBasis::Modal);
  const mesh::sf::SFTag basis_tag2b(ElemShape::Triag, SFunc::Lagrange, P4, ModalBasis::Modal);

  interpolation::PolyCacheDescription pdesc;
  const mesh::DiscreteElemKey key1 =
      mesh::DiscreteElemKey(elem_tag1_ext, basis_tag1, quad_tag1_ext);
  const mesh::DiscreteElemKey key2a =
      mesh::DiscreteElemKey(elem_tag2_ext, basis_tag2a, quad_tag2_ext);
  const mesh::DiscreteElemKey key2b =
      mesh::DiscreteElemKey(elem_tag2_ext, basis_tag2b, quad_tag2_ext);
  pdesc.add_to_block_count(key1, 4);
  pdesc.add_to_block_count(key2a, 4);
  pdesc.add_to_block_count(key2b, 2);

  pdesc.print();

  interpolation::PolyCache poly_cache;
  poly_cache.allocate(pdesc);
  BOOST_CHECK_EQUAL(poly_cache.nb_values_in_buffer(), 0u);
  BOOST_CHECK_EQUAL(poly_cache.capacity(), 10u);

  /*
  std::cout << "POLYNOMIAL CACHE CELLWISE VALUES:" << std::endl;
  for (Uint i = 0; i < 2; ++i)
  {
    poly_cache.push_back_to_buffer(key1);
    poly_cache.push_back_to_buffer(key2a);
  }

  for (Uint i = 0; i < 2; ++i)
  {
    poly_cache.push_back_to_buffer(key2b);
  }

  for (Uint i = 0; i < poly_cache.nb_values_in_buffer(); ++i)
  {
    const auto cell_vals = poly_cache.cellwise_values(i);
    std::cout << " ==>" << cell_vals.discrete_elem_key() << std::endl;
    std::cout << cell_vals.reference_sf_values() << std::endl;
  }

  for (Uint i = 0; i < poly_cache.nb_values_in_buffer(); ++i)
  {
    const auto cell_vals = poly_cache.cellwise_values(i);
    std::cout << " ==>" << cell_vals.discrete_elem_key() << std::endl;
    std::cout << "[Idx, pos_in_block] = [" << i << "," << poly_cache.position_in_block(i) << "]"
              << "\n"
              << std::endl;
  }
  */

  // BOOST_CHECK_EQUAL(poly_cache.nb_values_in_buffer(), 10u);
  // BOOST_CHECK_EQUAL(poly_cache.capacity(), 10u);

  // poly_cache.print();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(poly_metric_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  const Real tolerance = 1.e-7;
  test_metric<MeshConfig2D, _2D>(*(mesh2d->dof_storage("geo_dofs")), PointSetID::Gauss, 2.0,
                                 tolerance);
}

// ----------------------------------------------------------------------------
