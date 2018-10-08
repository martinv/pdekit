/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE geometry_metric_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "interpolation/FunctionSpace.hpp"
#include "interpolation/GeometryMetric.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------
// This function is supposed to be used with quadratures that place at
// least some of their points in element interior, such as standard
// Gauss quadrature
// ----------------------------------------------------------------------------
template <typename MeshConfig, Uint DIM>
void test_metric_inside_elems(
    typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
    typename interpolation::FunctionSpace<MeshConfig, DIM>::ptr fs_geometry,
    interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> &geo_metric,
    const PointSetID quadrature_type, Real reference_volume_value, const Real tol)
{
  clock_t t1, t2;
  Real elapsed;

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;

  std::set<mesh::PointSetTag> elem_types;
  for (Uint cg = 0; cg < cell_dofs.all_active_dof_groups().size(); ++cg)
  {
    const mesh::MeshEntity cell = (*cell_dofs.all_active_dof_groups()[cg].begin()).mesh_entity();
    elem_types.insert(cell.pt_set_id());
  }

  const auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  const auto eval_pt_set_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::PointSetTag(shape, P4, quadrature_type);
  };

  fs_geometry->set_reference_fe_values(
      common::make_iter_range(cell_dofs.cbegin(), cell_dofs.cend()), sf_generator,
      eval_pt_set_generator);

  geometry_cache.allocate(fs_geometry->discrete_elements().cbegin(),
                          fs_geometry->discrete_elements().cend(), cell_dofs.nb_active_cells());
  geo_metric.allocate_buffer(fs_geometry->discrete_elements().cbegin(),
                             fs_geometry->discrete_elements().cend(), cell_dofs.nb_active_cells());

  BOOST_CHECK_EQUAL(geo_metric.max_nb_blocks_in_buffer(), cell_dofs.nb_active_cells());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());
    const ElemShape cell_shape      = cell.pt_set_id().elem_shape();
    const Uint cell_order           = cell.pt_set_id().poly_order();
    const mesh::sf::SFTag basis_tag = sf_generator(cell_shape, cell_order);

    const mesh::PointSetTagExt support_pt_set_ext =
        mesh::PointSetTagExt(cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::PointSetTag eval_pt_set = eval_pt_set_generator(cell_shape, cell_order);
    const mesh::PointSetTagExt eval_pt_set_ext(eval_pt_set, P0, mesh::CellTransform::NO_TRANS, 0);

    const mesh::DiscreteElemKey geo_key(support_pt_set_ext, basis_tag, eval_pt_set_ext);

    geometry_cache.push_back_to_buffer(cell_coords, geo_key);
  }

  t1 = clock();

  geo_metric.evaluate(geometry_cache, interpolation::RebuildMetricIndex{true});

  t2 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (interpolation + computation of jacobians) = " << elapsed << " s"
            << std::endl;

  Real domain_volume = 0.0;

  typedef typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 DIM>::cellwise_metric cellwise_metric_type;

  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    cellwise_metric_type cell_met = geo_metric.cellwise_values(c);

    math::DenseConstVecView<Real> jdet   = cell_met.jdet();
    math::DenseDVec<Real> const &weights = cell_met.pt_weights();

    for (Uint q = 0; q < cell_met.nb_qd_pts(); ++q)
    {
      domain_volume += jdet[q] * weights[q];
    }
  }

  std::cout.precision(10);
  std::cout.setf(std::ios::fixed);
  std::cout << "Domain measure = " << std::setw(10) << domain_volume << " [expected "
            << reference_volume_value << "]" << std::endl;

  BOOST_CHECK_CLOSE(domain_volume, reference_volume_value, tol);
}

// ----------------------------------------------------------------------------
// This function is supposed to be used with quadratures that place points
// on element facets, such as FaceGauss quadrature
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void test_metric_on_elem_bdry(
    typename result_of::dof_map_t<MeshConfig> const &cell_dofs,
    typename interpolation::FunctionSpace<MeshConfig, DIM>::ptr fs_geometry,
    interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> &geo_metric,
    const PointSetID quadrature_type, Real reference_volume_value, const Real tol)
{
  clock_t t1, t2;
  Real elapsed;

  interpolation::GeometryCache<MeshConfig::GDIM> geometry_cache;

  const Uint cache_size = 10u;

  fs_geometry->set_reference_fe_values(cell_dofs, SFunc::Lagrange, P4, quadrature_type);
  geometry_cache.allocate(fs_geometry->reference_elements(), cache_size);

  std::cout << "TYPES ALLOWED IN GEOMETRY CACHE:" << std::endl;
  geometry_cache.print_types();

  geo_metric.allocate_buffer(fs_geometry->reference_elements(), cache_size);

  BOOST_CHECK_EQUAL(geo_metric.max_nb_blocks_in_buffer(), cache_size);

  t1 = clock();
  using cellwise_metric_type =
      typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                             DIM>::cellwise_metric;

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint cell_idx = 0; cell_idx < cell_dofs.nb_active_cells(); ++cell_idx)
  {
    geometry_cache.flush();
    geo_metric.empty_buffer();

    const mesh::CellTopologyView<MeshConfig> tcell_view = cell_dofs.tcell(cell_idx);
    const mesh::MeshEntity cell                         = cell_dofs.active_cell(cell_idx);

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());
    for (Uint sub_ent = 0; sub_ent < cell.nb_sub_elements(MeshConfig::TDIM - 1); ++sub_ent)
    {
      geometry_cache.push_back_to_buffer(
          cell_coords,
          mesh::PointSetTagExt(cell.pt_set_id(), P0, mesh::CellTransform::NO_TRANS, sub_ent));
    }

    geo_metric.evaluate(geometry_cache);

    for (Uint i = 0; i < geo_metric.nb_values_in_buffer(); ++i)
    {
      cellwise_metric_type cell_met = geo_metric.cellwise_values(i);

      /*
      std::cout << "[i = " << i << "]" << std::endl;
      std::cout << std::endl;
      std::cout << cell_met.coord_transf_derivatives(X0) << std::endl;
      std::cout << cell_met.coord_transf_derivatives(X1) << std::endl <<
      std::endl; for (Uint q = 0; q < cell_met.nb_qd_pts(); ++q)
      {
        std::cout << "\t" << cell_met.jdet()[q] << std::endl;
      }
      std::cout << std::endl;
      */
    }

  } // Loop over cells in mesh

  t2 = clock();

  elapsed = ((Real)(t2 - t1)) / CLOCKS_PER_SEC;
  std::cout << "CPU time (cellwise interpolation with face quadrature) = " << elapsed << " s"
            << std::endl;
}

// ----------------------------------------------------------------------------
// This function verifies that normals computed py Piola transform
// on element facets with FaceGauss quadrature are identical to normals
// computed by metric which uses purely facets
// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint DIM>
void test_metric_normals(mesh::Tria<MeshConfig> const &input_mesh,
                         typename result_of::dof_map_t<MeshConfig> const &cell_dofs)
{
  std::vector<mesh::DiscreteElemKey> ref_elems_restricted;
  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_restr_elems;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM> geo_metric_restr_elems;

  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_facets;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, DIM> geo_metric_facets;

  // FILL TYPE MAP WHICH USES ELEMENTS RESTRICTION TO SUB-ENTITY

  for (Uint f = 0; f < input_mesh.active_skeleton_size(DIM); ++f)
  {
    const mesh::TraceIncidences facet_block =
        input_mesh.active_skeleton_entry(DIM, mesh::ActiveIdx(f));

    for (Uint cell_idx_in_block = 0; cell_idx_in_block < facet_block.size(); ++cell_idx_in_block)
    {
      const mesh::CellTopologyView<MeshConfig> tcell =
          input_mesh.cell(mesh::FlatIdx(facet_block.cell_id(cell_idx_in_block)));
      const mesh::MeshEntity elem_entity =
          cell_dofs.active_cell(mesh::ActiveIdx(tcell.active_idx()));
      const mesh::PointSetTag std_reg_tag_elem(elem_entity.pt_set_id());

      const Uint sub_ent = facet_block.local_id(cell_idx_in_block);

      const mesh::PointSetTagExt std_reg_tag_elem_ext(
          std_reg_tag_elem, P1, mesh::CellTransform::RESTRICT_TO_CODIM_1, sub_ent);

      const mesh::sf::SFTag sf_type(std_reg_tag_elem.elem_shape(), SFunc::Lagrange,
                                    std_reg_tag_elem.poly_order(), ModalBasis::Modal);

      mesh::PointSetTag quad_tag(std_reg_tag_elem.elem_shape(), P4, PointSetID::FaceGauss);
      mesh::PointSetTagExt quad_tag_ext(quad_tag, P1, mesh::CellTransform::RESTRICT_TO_CODIM_1,
                                        sub_ent);

      const mesh::DiscreteElemKey key(std_reg_tag_elem_ext, sf_type, quad_tag_ext);

      bool key_is_registered = false;
      for (const auto &elem : ref_elems_restricted)
      {
        if (elem == key)
        {
          key_is_registered = true;
        }
      }
      if (!key_is_registered)
      {
        ref_elems_restricted.push_back(key);
      }
    }
  } // Loop over mesh facets

  // FILL TYPE MAP WHICH USES DIRECTLY ENTITIES OF CODIMENSION 1 (ELEMENT
  // TRACE ENTITIES)

  std::vector<mesh::DiscreteElemKey> ref_elems_on_traces;

  for (Uint f = 0; f < input_mesh.active_skeleton_size(DIM); ++f)
  {
    const mesh::TraceIncidences facet_block =
        input_mesh.active_skeleton_entry(DIM, mesh::ActiveIdx(f));

    for (Uint cell_idx_in_block = 0; cell_idx_in_block < facet_block.size(); ++cell_idx_in_block)
    {
      const mesh::CellTopologyView<MeshConfig> tcell =
          input_mesh.cell(mesh::FlatIdx(facet_block.cell_id(cell_idx_in_block)));
      const mesh::MeshEntity elem_entity =
          cell_dofs.active_cell(mesh::ActiveIdx(tcell.active_idx()));

      const Uint sub_ent = facet_block.local_id(cell_idx_in_block);

      const mesh::MeshEntity facet_entity = elem_entity.sub_entity(DIM, sub_ent);

      const mesh::PointSetTag std_reg_tag_facet(facet_entity.pt_set_id());

      const mesh::PointSetTagExt std_reg_tag_facet_ext(std_reg_tag_facet, P1,
                                                       mesh::CellTransform::NO_TRANS, 0);

      const mesh::sf::SFTag sf_type(std_reg_tag_facet.elem_shape(), SFunc::Lagrange,
                                    std_reg_tag_facet.poly_order(), ModalBasis::Modal);

      const mesh::PointSetTag quad_tag(std_reg_tag_facet.elem_shape(), P4, PointSetID::Gauss);
      const mesh::PointSetTagExt quad_tag_ext(quad_tag, P1, mesh::CellTransform::NO_TRANS, 0);

      const mesh::DiscreteElemKey key(std_reg_tag_facet_ext, sf_type, quad_tag_ext);

      bool key_is_registered = false;

      for (const auto &elem : ref_elems_on_traces)
      {
        if (elem == key)
        {
          key_is_registered = true;
        }
      }

      if (!key_is_registered)
      {
        ref_elems_on_traces.push_back(key);
      }
    }
  } // Loop over mesh facets

  geo_cache_restr_elems.allocate(ref_elems_restricted.cbegin(), ref_elems_restricted.cend(), 1);
  geo_metric_restr_elems.allocate_buffer(ref_elems_restricted.cbegin(), ref_elems_restricted.cend(),
                                         1);

  geo_cache_facets.allocate(ref_elems_on_traces.cbegin(), ref_elems_on_traces.cend(), 1);
  geo_metric_facets.allocate_buffer(ref_elems_on_traces.cbegin(), ref_elems_on_traces.cend(), 1);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint f = 0; f < input_mesh.active_skeleton_size(DIM); ++f)
  {
    const mesh::TraceIncidences facet_block =
        input_mesh.active_skeleton_entry(DIM, mesh::ActiveIdx(f));

    for (Uint cell_idx_in_block = 0; cell_idx_in_block < facet_block.size(); ++cell_idx_in_block)
    {
      const mesh::CellTopologyView<MeshConfig> tcell =
          input_mesh.cell(mesh::FlatIdx(facet_block.cell_id(cell_idx_in_block)));

      const mesh::MeshEntity elem_entity =
          cell_dofs.active_cell(mesh::ActiveIdx(tcell.active_idx()));

      const Uint sub_ent = facet_block.local_id(cell_idx_in_block);

      const mesh::MeshEntity facet_entity = elem_entity.sub_entity(DIM, sub_ent);

      const mesh::PointSetTag std_reg_tag_elem(elem_entity.pt_set_id());
      const mesh::PointSetTagExt std_reg_tag_elem_ext(
          std_reg_tag_elem, P1, mesh::CellTransform::RESTRICT_TO_CODIM_1, sub_ent);

      const mesh::sf::SFTag sf_type_elem(std_reg_tag_elem.elem_shape(), SFunc::Lagrange,
                                         std_reg_tag_elem.poly_order(), ModalBasis::Modal);

      mesh::PointSetTag quad_tag_elem(std_reg_tag_elem.elem_shape(), P4, PointSetID::FaceGauss);
      mesh::PointSetTagExt quad_tag_elem_ext(quad_tag_elem, P1,
                                             mesh::CellTransform::RESTRICT_TO_CODIM_1, sub_ent);

      const mesh::DiscreteElemKey key_elem(std_reg_tag_elem_ext, sf_type_elem, quad_tag_elem_ext);

      const mesh::PointSetTag std_reg_tag_facet(facet_entity.pt_set_id());
      const mesh::PointSetTagExt std_reg_tag_facet_ext(std_reg_tag_facet, P1,
                                                       mesh::CellTransform::NO_TRANS, 0);
      const mesh::sf::SFTag sf_type_facet(std_reg_tag_facet.elem_shape(), SFunc::Lagrange,
                                          std_reg_tag_facet.poly_order(), ModalBasis::Modal);

      const mesh::PointSetTag quad_tag_facet(std_reg_tag_facet.elem_shape(), P4, PointSetID::Gauss);
      const mesh::PointSetTagExt quad_tag_facet_ext(quad_tag_facet, P1,
                                                    mesh::CellTransform::NO_TRANS, 0);

      const mesh::DiscreteElemKey key_facet(std_reg_tag_facet_ext, sf_type_facet,
                                            quad_tag_facet_ext);

      geo_cache_restr_elems.flush();
      geo_cache_facets.flush();

      const math::DenseConstMatView<Real> elem_coords = loc_interpolator.transfer_coords(
          tcell.pt_set_id(), elem_entity.pt_set_id(), tcell.coordinates());

      geo_cache_restr_elems.push_back_to_buffer(elem_coords, key_elem);

      const math::DenseConstMatView<Real> facet_coords = loc_interpolator.transfer_coords(
          tcell.pt_set_id(DIM, sub_ent), facet_entity.pt_set_id(), tcell.coordinates(DIM, sub_ent));
      geo_cache_facets.push_back_to_buffer(facet_coords, key_facet);

      geo_metric_restr_elems.empty_buffer();
      geo_metric_restr_elems.evaluate(geo_cache_restr_elems,
                                      interpolation::RebuildMetricIndex{true});

      geo_metric_facets.empty_buffer();
      geo_metric_facets.evaluate(geo_cache_facets, interpolation::RebuildMetricIndex{true});

      const typename interpolation::GeometryMetric<MeshConfig::GDIM,
                                                   MeshConfig::TDIM>::cellwise_metric geo_cell_met =
          geo_metric_restr_elems.cellwise_values(0);

      const typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                   DIM>::cellwise_metric geo_facet_met =
          geo_metric_facets.cellwise_values(0);

      /*
      std::cout << "**************************************************" <<
      std::endl;

      std::cout << elem_entity << std::endl;
      std::cout << facet_entity << std::endl;

      std::cout << elem_coords << std::endl;
      std::cout << "Local entity = " << sub_ent << std::endl;
      std::cout << facet_coords << std::endl;

      std::cout << "Interpolated coords on element: " << std::endl;
      std::cout << geo_cell_met.interpolated_coords() << std::endl;

      std::cout << "Interpolated coords on facet: " << std::endl;
      std::cout << geo_facet_met.interpolated_coords() << std::endl;


      std::cout << "Normals on element: " << std::endl;
      std::cout << geo_cell_met.normals() << std::endl;

      std::cout << "Normals on facet: " << std::endl;
      std::cout << geo_facet_met.normals() << std::endl;

      std::cout << "**************************************************" <<
      std::endl;
      */

      math::DenseConstMatView<Real> normals_cell_restr = geo_cell_met.normals();
      math::DenseConstMatView<Real> normals_facet      = geo_facet_met.normals();

      BOOST_CHECK_EQUAL(normals_cell_restr.rows(), normals_facet.rows());

      for (Uint r = 0; r < normals_cell_restr.rows(); ++r)
      {
        Real max_norm = 0.0;
        for (Uint d = 0; d < MeshConfig::GDIM; ++d)
        {
          max_norm += std::abs(normals_cell_restr(r, d) - normals_facet(r, d));
        }

        BOOST_CHECK_LT(max_norm, 1.e-10);
      }
    }
  } // Loop over mesh facets
}

// ----------------------------------------------------------------------------

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", *mesh2d, "geo_dofs");

  interpolation::FunctionSpace<MeshConfig2D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D>>();

  const Real tolerance = 1.e-7;
  interpolation::GeometryMetric<MeshConfig2D::GDIM, MeshConfig2D::TDIM, _2D> geo_metric;
  test_metric_inside_elems<MeshConfig2D, _2D>(*(mesh2d->dof_storage("geo_dofs")), fs_geometry,
                                              geo_metric, PointSetID::Gauss, 1.0, tolerance);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_3D_utest)
{
  MeshType3D::shared_ptr mesh3d = std::make_shared<MeshType3D>("mesh3D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  interpolation::FunctionSpace<MeshConfig3D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig3D>>();

  const Real tolerance = 1.e-7;
  interpolation::GeometryMetric<MeshConfig3D::GDIM, MeshConfig3D::TDIM, _3D> geo_metric;
  test_metric_inside_elems<MeshConfig3D, _3D>(*(mesh3d->dof_storage("geo_dofs")), fs_geometry,
                                              geo_metric, PointSetID::Gauss, 1.0, tolerance);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_normals_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  // mesh_reader.read_mesh_from_file("unit_square_mixed_p2.msh", *mesh2d,
  // "geo_dofs");
  mesh_reader.read_mesh_from_file("unit_circle_tri_p3.msh", *mesh2d, "geo_dofs");

  test_metric_normals<MeshConfig2D, _1D>(*mesh2d, *(mesh2d->dof_storage("geo_dofs")));
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_normals_3D_utest)
{
  MeshType3D::shared_ptr mesh3d = std::make_shared<MeshType3D>("mesh2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  test_metric_normals<MeshConfig3D, _2D>(*mesh3d, *(mesh3d->dof_storage("geo_dofs")));
}

// ----------------------------------------------------------------------------

#if 0
BOOST_AUTO_TEST_CASE(metric_1D_boundary_utest)
{
  MeshType2D::shared_ptr mesh2d = MeshType2D::shared_ptr(new MeshType2D("circle2D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_circle_tri_p3.msh", *mesh2d);

  mesh::MeshBoundarySet<MeshConfig2D> const &boundary = mesh2d->topology().all_boundaries();

  interpolation::FunctionSpace<MeshConfig2D, _1D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D, _1D>>(mesh2d);

  const Real tolerance = 1.e-5;
  interpolation::GeometryMetric<MeshConfig2D, _1D> geo_metric;
  test_metric_inside_elems<MeshConfig2D, _1D>(*boundary.domain("boundary"), mesh2d->geometry(),
                                              fs_geometry, geo_metric, Gauss, 2. * PI, tolerance);
}
#endif

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(metric_2D_mesh_1D_facets_utest)
{
#if 0
  MeshType2D::shared_ptr mesh2d = MeshType2D::shared_ptr(new MeshType2D("square_regular2D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2_regular.msh", *mesh2d);

  mesh::TopologyAlgorithms::compute_connectivity(mesh2d->topology(), _1D, _0D, true);

  const result_of::cells<MeshConfig2D>::type &facets = mesh2d->topology().cells(_1D);

  typedef result_of::geometry<MeshConfig2D>::type geometry_type;
  const geometry_type &mesh_geometry = mesh2d->geometry();

  interpolation::FunctionSpace<MeshConfig2D, _1D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D, _1D>>(mesh2d);

  const Real tolerance = 1.e-12;

  // The unit square domain is meshed with 16 x 16 small squares, each of which
  // is divided by a diagonal into 2 triangles. The total length of edges in
  // the mesh is
  // 17 (total length of all horizontal edges)
  // +
  // 17 (total length of all vertical edges)
  // +
  //      1/16 * sqrt(2)      x               (16 * 16)
  // (length of one diagonal) x (number of squares, i.e. triangle pairs)

  const Real total_edge_length = 17. + 17. + (1. / 16. * std::sqrt(2)) * (16. * 16.);
  interpolation::GeometryMetric<MeshConfig2D, _1D> geo_metric;
  test_metric_inside_elems<MeshConfig2D, _1D>(facets, mesh_geometry, fs_geometry, geo_metric, Gauss,
                                              total_edge_length, tolerance);

  // Write vectors to file and plot them with gnuplot as
  //  plot ’file.dat’ using 1:2:3:4 with vectors head filled lt 2
  /*
      plot 'facets_2D.dat' using 1:2:3:4 with vectors, \
           'facet_normals_2D.dat' using 1:2:3:4 with vectors
  */
  typedef interpolation::GeometryMetric<MeshConfig2D, _1D>::cellwise_metric cellwise_metric_type;

  std::ofstream outfile;
  outfile.setf(std::ios::fixed);
  outfile.precision(10);

  outfile.open("facets_2D.dat");

  math::StaticVector<Real, 2> diff;

  for (Uint c = 0; c < facets.nb_cells(); ++c)
  {
    const mesh::MeshEntity facet = facets.cell(c);
    const geometry_type::point_coord pt_start = mesh_geometry.node(facet.vertex(0));
    const geometry_type::point_coord pt_end = mesh_geometry.node(facet.vertex(1));

    diff = pt_end - pt_start;
    outfile << pt_start[X0] << " " << pt_start[X1] << " " << diff << std::endl;
  }

  outfile.close();

  outfile.open("facet_normals_2D.dat");
  for (Uint c = 0; c < facets.nb_cells(); ++c)
  {
    cellwise_metric_type cell_met = geo_metric.cellwise_values(c);

    math::ConstMatrixBlock<Real> quad_pt_coords = cell_met.interpolated_coords();
    math::ConstMatrixBlock<Real> quad_pt_normals = cell_met.normals();

    for (Uint pt = 0; pt < cell_met.nb_qd_pts(); ++pt)
    {
      const math::ConstVectorBlock<Real> normal = quad_pt_normals.row_transpose(pt);
      // The factor 1/32 is here just to scale the normals for plotting
      outfile << quad_pt_coords.row(pt) << " " << 1. / 32. * normal[X0] << " "
              << 1. / 32. * normal[X1] << std::endl;
    }
  }
  outfile.close();
#endif
}

// ----------------------------------------------------------------------------

#if 0
BOOST_AUTO_TEST_CASE(metric_2D_boundary_utest)
{
  MeshType3D::shared_ptr mesh3d = MeshType3D::shared_ptr(new MeshType3D("sphere2d"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_sphere_tet_p3.msh", *mesh3d);

  mesh::MeshBoundarySet<MeshConfig3D> const &boundary = mesh3d->topology().all_boundaries();

  interpolation::FunctionSpace<MeshConfig3D, _2D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig3D, _2D>>(mesh3d);

  const Real tolerance = 7.e-5;
  interpolation::GeometryMetric<MeshConfig3D, _2D> geo_metric;
  test_metric_inside_elems<MeshConfig3D, _2D>(*boundary.domain("boundary"), mesh3d->geometry(),
                                              fs_geometry, geo_metric, Gauss, 4. * PI, tolerance);
}
#endif

// ----------------------------------------------------------------------------

/*
BOOST_AUTO_TEST_CASE(metric_2D_mesh_face_quadrature_utest)
{
  MeshType2D::shared_ptr mesh2d = MeshType2D::shared_ptr(new
MeshType2D("circle2D"));

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_sphere_tet_p2.msh", *mesh2d,
"geo_dofs");

  interpolation::FunctionSpace<MeshConfig2D, _2D>::ptr fs_geometry =
      std::make_shared<interpolation::FunctionSpace<MeshConfig2D, _2D>>();

  const Real tolerance = 7.e-5;
  interpolation::GeometryMetric<MeshConfig2D, _2D> geo_metric;
  test_metric_on_elem_bdry<MeshConfig2D,
_2D>(*(mesh2d->dof_storage("geo_dofs")), fs_geometry, geo_metric,
PointSetID::FaceGauss, 4. * PI, tolerance);
}
*/

// ----------------------------------------------------------------------------
