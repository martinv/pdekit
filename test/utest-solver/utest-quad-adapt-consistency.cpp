/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_adapt_consistency_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"
#include "tools/MeshInspector.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void check_adapted_mesh(const Tria<MeshConfig> &mesh, const std::string &dof_handler_name,
                        const Uint quad_order)
{
  const typename result_of::dof_map_t<MeshConfig> &cell_dofs =
      *(mesh.dof_storage(dof_handler_name));

  std::vector<Real> cell_volume(cell_dofs.nb_active_cells());
  cell_volume.assign(cell_dofs.nb_active_cells(), 0.0);

  interpolation::FunctionSpace<MeshConfig, MeshConfig::TDIM - 1> skeleton_space;

  auto sf_generator = [=](const ElemShape shape, const Uint order) {
    return mesh::sf::SFTag(shape, SFunc::Lagrange, order, ModalBasis::Modal);
  };

  auto quad_generator = [quad_order](const ElemShape shape, const Uint elem_order) {
    return mesh::PointSetTag(shape, quad_order, PointSetID::Gauss);
  };

  skeleton_space.set_reference_fe_values_on_skeleton(mesh, cell_dofs, sf_generator, quad_generator);
  // skeleton_space.print_element_types();

  mesh::QuadraturePermutation quad_p_L, quad_p_R;

  /*
  for (Uint c = 0; c < cell_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = cell_dofs.active_cell(c);
    std::cout << "cell [" << cell.idx() << "] = " << cell << std::endl;
  }
  */

  const Uint facet_dim = MeshConfig::TDIM - 1;

  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_L;
  interpolation::GeometryCache<MeshConfig::GDIM> geo_cache_R;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, facet_dim> geo_metric_L;
  interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM, facet_dim> geo_metric_R;

  geo_cache_L.allocate(skeleton_space.discrete_elements().cbegin(),
                       skeleton_space.discrete_elements().cend(), 1);
  geo_cache_R.allocate(skeleton_space.discrete_elements().cbegin(),
                       skeleton_space.discrete_elements().cend(), 1);

  geo_metric_L.allocate_buffer(skeleton_space.discrete_elements().cbegin(),
                               skeleton_space.discrete_elements().cend(), 1);
  geo_metric_R.allocate_buffer(skeleton_space.discrete_elements().cbegin(),
                               skeleton_space.discrete_elements().cend(), 1);

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint f = 0; f < mesh.active_skeleton_size(); ++f)
  {
    const mesh::TraceIncidences facet_block =
        mesh.active_skeleton_entry(facet_dim, mesh::ActiveIdx(f));

    const Uint local_idx_L = 0;
    const Uint local_idx_R = (facet_block.size() == 2) ? 1 : 0;

    // NOTE THAT FACET BLOCK ONLY KNOWS __ABSOLUTE__ (LINEAR) POSITIONS OF
    // CELLS, NOT THEIR ACTIVE INDICES !!!
    const mesh::CellTopologyView<MeshConfig> tcell_L =
        mesh.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_L)));
    const mesh::CellTopologyView<MeshConfig> tcell_R =
        mesh.cell(mesh::FlatIdx(facet_block.cell_id(local_idx_R)));

    /*
    std::cout <<
    "----------------------------------------------------------------" <<
    std::endl; std::cout << "Topology cell left: active = " <<
    tcell_L.active_idx()
              << ", absolute = " << tcell_L.linear_pos_idx() << std::endl;
    std::cout << "Block entity left: active cell id = "
              << facet_block.active_cell_id(local_idx_L) << std::endl;
    std::cout << "Block entity left: cell id = " <<
    facet_block.cell_id(local_idx_L) << std::endl
              << std::endl;
    std::cout << "Topology cell right: active = " << tcell_R.active_idx()
              << ", absolute = " << tcell_R.linear_pos_idx() << std::endl;
    std::cout << "Block entity right: active cell id = "
              << facet_block.active_cell_id(local_idx_R) << std::endl;
    std::cout << "Block entity right: cell id = " <<
    facet_block.cell_id(local_idx_R)
              << std::endl
              << std::endl;
    */

    const ActiveIdx active_cell_id_L = tcell_L.active_idx();
    const ActiveIdx active_cell_id_R = tcell_R.active_idx();

    mesh::MeshEntity facet_L = cell_dofs.active_cell(active_cell_id_L);
    facet_L.local_transform(facet_dim, facet_block.local_id(local_idx_L));

    mesh::MeshEntity facet_R = cell_dofs.active_cell(active_cell_id_R);
    facet_R.local_transform(facet_dim, facet_block.local_id(local_idx_R));

    const mesh::EntityRealignCode local_sub_tag_L =
        facet_block.permutation(local_idx_L).get().code();

    const ElemShape facet_shape_L = facet_L.pt_set_id().elem_shape();
    const Uint facet_order_L      = facet_L.pt_set_id().poly_order();
    const mesh::PointSetTagExt pt_set_facet_ext_L(facet_L.pt_set_id(), P0,
                                                  local_sub_tag_L.adapt_op_id(),
                                                  local_sub_tag_L.local_pos_in_parent());
    const mesh::sf::SFTag sf_facet_L     = sf_generator(facet_shape_L, facet_order_L);
    const mesh::PointSetTag quad_facet_L = quad_generator(facet_shape_L, facet_order_L);
    const mesh::PointSetTagExt quad_facet_ext_L(quad_facet_L, P0, local_sub_tag_L.adapt_op_id(),
                                                local_sub_tag_L.local_pos_in_parent());

    const mesh::DiscreteElemKey key_L(pt_set_facet_ext_L, sf_facet_L, quad_facet_ext_L);

    const mesh::EntityRealignCode local_sub_tag_R =
        facet_block.permutation(local_idx_R).get().code();

    const ElemShape facet_shape_R = facet_R.pt_set_id().elem_shape();
    const Uint facet_order_R      = facet_R.pt_set_id().poly_order();
    const mesh::PointSetTagExt pt_set_facet_ext_R(facet_R.pt_set_id(), P0,
                                                  local_sub_tag_R.adapt_op_id(),
                                                  local_sub_tag_R.local_pos_in_parent());
    const mesh::sf::SFTag sf_facet_R     = sf_generator(facet_shape_R, facet_order_R);
    const mesh::PointSetTag quad_facet_R = quad_generator(facet_shape_R, facet_order_R);
    const mesh::PointSetTagExt quad_facet_ext_R(quad_facet_R, P0, local_sub_tag_R.adapt_op_id(),
                                                local_sub_tag_R.local_pos_in_parent());

    const mesh::DiscreteElemKey key_R(pt_set_facet_ext_R, sf_facet_R, quad_facet_ext_R);

    geo_cache_L.flush();
    geo_cache_R.flush();

    const PointSetTag tcell_facet_L =
        tcell_L.pt_set_id(facet_dim, facet_block.local_id(local_idx_L));

    const math::DenseConstMatView<Real> geo_facet_coords_L = loc_interpolator.transfer_coords(
        tcell_facet_L, facet_L.pt_set_id(),
        tcell_L.coordinates(facet_dim, facet_block.local_id(local_idx_L)));

    geo_cache_L.push_back_to_buffer(geo_facet_coords_L, key_L);
    geo_metric_L.empty_buffer();
    geo_metric_L.evaluate(geo_cache_L, interpolation::RebuildMetricIndex{true});

    const PointSetTag tcell_facet_R =
        tcell_R.pt_set_id(facet_dim, facet_block.local_id(local_idx_R));

    const math::DenseConstMatView<Real> geo_facet_coords_R = loc_interpolator.transfer_coords(
        tcell_facet_R, facet_R.pt_set_id(),
        tcell_R.coordinates(facet_dim, facet_block.local_id(local_idx_R)));

    geo_cache_R.push_back_to_buffer(geo_facet_coords_R, key_R);
    geo_metric_R.empty_buffer();
    geo_metric_R.evaluate(geo_cache_R, interpolation::RebuildMetricIndex{true});

    // std::cout << key_L << std::endl;
    const typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 facet_dim>::cellwise_metric facet_met_data_L =
        geo_metric_L.cellwise_values(0);
    const math::DenseConstMatView<Real> qcoords_L = facet_met_data_L.interpolated_coords();
    // std::cout << qcoords_L << std::endl;
    const math::DenseConstVecView<Real> jdet_L = facet_met_data_L.jdet();
    const math::DenseDVec<Real> &wq_L          = facet_met_data_L.pt_weights();

    // std::cout << facet_met_data_L.interpolated_coords() << std::endl;

    // std::cout << key_R << std::endl;
    const typename interpolation::GeometryMetric<MeshConfig::GDIM, MeshConfig::TDIM,
                                                 facet_dim>::cellwise_metric facet_met_data_R =
        geo_metric_R.cellwise_values(0);
    const math::DenseConstMatView<Real> qcoords_R = facet_met_data_R.interpolated_coords();
    // std::cout << qcoords_R << std::endl;
    const math::DenseConstVecView<Real> jdet_R = facet_met_data_R.jdet();
    const math::DenseDVec<Real> &wq_R          = facet_met_data_R.pt_weights();

    BOOST_CHECK_EQUAL(qcoords_L.rows(), qcoords_R.rows());

    const mesh::PointSetTag contour_quad_tag_L =
        mesh::PointSetTag(facet_L.pt_set_id().elem_shape(), quad_order, PointSetID::Gauss);

    const mesh::PointSetTag contour_quad_tag_R =
        mesh::PointSetTag(facet_R.pt_set_id().elem_shape(), quad_order, PointSetID::Gauss);

    quad_p_L.change_type(contour_quad_tag_L, facet_block.local_id(local_idx_L),
                         facet_block.permutation(local_idx_L).get().code());

    quad_p_R.change_type(contour_quad_tag_R, facet_block.local_id(local_idx_R),
                         facet_block.permutation(local_idx_R).get().code());

    math::DenseConstMatView<Real> const normals_L = facet_met_data_L.normals();
    math::DenseConstMatView<Real> const normals_R = facet_met_data_R.normals();

    // std::cout << " ================================================ " <<
    // std::endl;

    /*
    std::cout << "Checking topological cells with absolute positions "
              << tcell_L.linear_pos_idx() << " and " <<
    tcell_R.linear_pos_idx()
              << std::endl;
    */

    Real length_L = 0;
    Real length_R = 0;

    for (Uint q = 0; q < qcoords_L.rows(); ++q)
    {
      // (Local) index of left quadrature point
      const Uint qV_L = quad_p_L.get().vertex(q);
      // (Local) index of right quadrature point
      const Uint qV_R = quad_p_R.get().vertex(q);

      /*
      std::cout << "q = " << q << std::endl;
      std::cout << "jL = " << jdet_L[qV_L] << std::endl;
      std::cout << "wL = " << wq_L[qV_L] << std::endl;
      std::cout << "nL = " << normals_L.row(qV_L) << std::endl;
      std::cout << "XL = " << qcoords_L.row(qV_L) << std::endl;

      std::cout << " <<<< " << std::endl;

      std::cout << "jR = " << jdet_R[qV_R] << std::endl;
      std::cout << "wR = " << wq_R[qV_R] << std::endl;
      std::cout << "nR = " << normals_R.row(qV_R) << std::endl;
      std::cout << "XR = " << qcoords_R.row(qV_R) << std::endl;
      */

      const Real one_over_dim = 1. / MeshConfig::GDIM;

      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {

        // Accumulate to the cell volume on the left:
        cell_volume[active_cell_id_L.id()] +=
            one_over_dim * jdet_L[qV_L] * wq_L[qV_L] * (normals_L(qV_L, d) * qcoords_L(qV_L, d));

        if (facet_block.size() == 2)
        {
          // Accumulate to the cell volume on the right:
          cell_volume[active_cell_id_R.id()] +=
              one_over_dim * jdet_R[qV_R] * wq_R[qV_R] * (normals_R(qV_R, d) * qcoords_R(qV_R, d));
        }

        BOOST_CHECK_CLOSE(qcoords_L(qV_L, d), qcoords_R(qV_R, d), 1.e-10);
      }
      length_L += jdet_L[qV_L] * wq_L[qV_L];

      if (facet_block.size() == 2)
      {
        length_R += jdet_R[qV_R] * wq_L[qV_R];
      }
    }

    /*
    std::cout << "   length L = " << length_L << std::endl;

    if (facet_block.size() == 2)
    {
      std::cout << "   length R = " << length_R << std::endl;
    }
    */
  }

  Real tot_volume = 0.0;
  for (Uint i = 0; i < cell_volume.size(); ++i)
  {
    tot_volume += cell_volume[i];
  }

  /*
  std::cout << "Cell volumes: " << std::endl;
  for (Uint i = 0; i < cell_volume.size(); ++i)
  {
    std::cout << "Active id[" << i << "]: volume = " << cell_volume[i] <<
  std::endl;
  }
  */

  std::cout << "Total evaluated domain volume = " << tot_volume << std::endl;
}

// ----------------------------------------------------------------------------

template <Uint GeoDim>
bool circle_intersects_element(const Real xc, const Real yc, const Real R,
                               const CellGeometry<GeoDim> &cell_coords)
{
  const Uint n_pts = cell_coords.size();

  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    const math::DenseVecView<const Real> point_coords0 = cell_coords.const_node_view(n);
    const math::DenseVecView<const Real> point_coords1 =
        cell_coords.const_node_view((n + 1) % n_pts);

    const Real x0 = point_coords0[X0];
    const Real y0 = point_coords0[X1];
    const Real x1 = point_coords1[X0];
    const Real y1 = point_coords1[X1];

    const Real dx = x1 - x0;
    const Real dy = y1 - y0;

    const Real A = dx * dx + dy * dy;
    const Real B = 2.0 * dx * (x0 - xc) + 2.0 * dy * (y0 - yc);
    const Real C = (x0 - xc) * (x0 - xc) + (y0 - yc) * (y0 - yc) - R * R;

    const Real det = B * B - 4.0 * A * C;

    if (det >= 0.0)
    {
      const Real t0 = (-B - std::sqrt(det)) / (2.0 * A);
      const Real t1 = (-B + std::sqrt(det)) / (2.0 * A);

      if ((0.0 <= t0) && (t0 <= 1.0))
        return true;

      if ((0.0 <= t1) && (t1 <= 1.0))
        return true;
    }
  }

  return false;
}

// ----------------------------------------------------------------------------

struct AdaptQuadratureConsistencyUtestFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(AdaptQuadConsistency_TestSuite, AdaptQuadratureConsistencyUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(h_adapt_select_cells_for_adaptation_w_hanging_nodes)
{
  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  const Uint nb_cells = mesh2d.nb_active_cells();

  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2d.adapt(adapt_schedule);

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));
  dof_handler.adapt(adapt_schedule);

  // Mark cells for new refinement
  cell_adapt_ops.resize(dof_handler.nb_active_cells());
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8]  = CellTransform::UNIFORM_REFINE;
  cell_adapt_ops[17] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

  interpolation::ScalarMeshFunction<Uint> marked_cells("", "marked_cells");
  marked_cells.resize(dof_handler.nb_active_cells());

  for (Uint i = 0; i < cell_adapt_ops.size(); ++i)
  {
    if (cell_adapt_ops[i] == CellTransform::UNIFORM_REFINE)
    {
      marked_cells[i] = 1;
    }
    else
    {
      marked_cells[i] = 0;
    }
  }

  // Write the mesh in the output file
  const std::string outfilename2d = "cells_marked_for_refinement_hanging_nodes.msh";
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, marked_cells, "adaptation marker",
                                          0, 0);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize(dof_handler.nb_active_cells());
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "active_cell_id");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(h_adapt_select_cells_for_adaptation_red_green)
{
  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  const Uint nb_cells = mesh2d.nb_active_cells();

#if 0
  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;
  // cell_adapt_ops[14] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);
  mesh2d.adapt(adapt_schedule);
  std::cout << ">>> TOPOLOGY ADAPTATION FINISHED" << std::endl;

  result_of::dof_handler<Cart2D>::type &dof_handler = *(mesh2d.dof_storage("geo_dofs"));
  dof_handler.adapt(adapt_schedule);

  std::cout << ">>> ADAPTATION FINISHED" << std::endl;

  std::vector<Uint> colors;

  // Re-set the adaptation operation vector
  cell_adapt_ops.resize(dof_handler.nb_active_cells());
  cell_adapt_ops.assign(dof_handler.nb_active_cells(), CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;

#else

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));
  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;
  // cell_adapt_ops[14] = CellTransform::UNIFORM_REFINE;

  std::vector<Uint> colors;
  colors.resize(dof_handler.nb_active_cells());
#endif

  mesh2d.color_cells_red_green(cell_adapt_ops, colors);

  interpolation::ScalarMeshFunction<Uint> marked_cells("", "marked_cells");
  marked_cells.resize(dof_handler.nb_active_cells());

  for (Uint i = 0; i < cell_adapt_ops.size(); ++i)
  {
    if (cell_adapt_ops[i] != CellTransform::NO_TRANS)
    {
      marked_cells[i] = 1;
    }
    else
    {
      marked_cells[i] = 0;
    }
  }

  interpolation::ScalarMeshFunction<Uint> cell_colors("", "cell_colors");
  cell_colors.resize(colors.size());
  for (Uint i = 0; i < colors.size(); ++i)
  {
    cell_colors[i] = colors[i];
  }

  std::cout << "Marked cells size = " << marked_cells.nb_entries() << std::endl;
  std::cout << "Cell colors size = " << cell_colors.nb_entries() << std::endl;
  std::cout << "Colors size = " << colors.size() << std::endl;

  // Write the mesh in the output file
  const std::string outfilename2d = "cells_marked_for_refinement_red_green.msh";
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, marked_cells, "adaptation marker",
                                          0, 0);

  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_colors, "cell_colors", 0, 0);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize(dof_handler.nb_active_cells());
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "active_cell_id");

  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "linear_cell_id");

  /*
  const CellTopologyView<Cart2D> center = mesh2d.active_cell(23);
  const std::vector<std::tuple<CellTopologyView<Cart2D>, Uint, Uint>> neighbours =
      mesh2d.active_neighbours(center);

  std::cout << "Neibhgours of active cell nr. " << center.active_idx() <<
  std::endl; for (auto neighb_item : neighbours)
  {
    const CellTopologyView<Cart2D> ncell = std::get<0>(neighb_item);
    const Uint n_loc_id = std::get<1>(neighb_item);
    const Uint center_loc_id = std::get<2>(neighb_item);

    std::cout << "Cell: [" << ncell.active_idx() << "," << n_loc_id << "]  ";
    std::cout << "local id in center cell: " << center_loc_id << std::endl;
  }
  */

  mesh2d.write_dual_graph_to_gmsh("geo_dofs", "dual.msh");

  // adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops,
  // h_AdaptStrategy::red_green); mesh2d.adapt(adapt_schedule);
  /*
  dof_handler.adapt(adapt_schedule);
  */
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(adapted_quadrature_consistency_utest)
{
  std::cout << "\n************ RUNNING MESH h-ADAPTATION TEST - 2D ************" << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  // const std::string infilename2d = "unit_square_mini.msh";
  const std::string infilename2d = "unit_square_mini_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  const std::vector<std::vector<Uint>> active_cells_to_refine = {{0}, {6}, {0, 1, 3, 5}};

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));

  MeshAdaptSequence<Cart2D> adapt_schedule;

  const Uint nb_adapt_passes = 6;
  std::vector<CellTransform> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    std::cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << std::endl;
    std::cout << "Adaptation pass #" << p + 1 << std::endl;
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    if (p < active_cells_to_refine.size())
    {
      for (Uint i = 0; i < active_cells_to_refine[p].size(); ++i)
      {
        cell_adapt_ops[active_cells_to_refine[p][i]] = CellTransform::UNIFORM_REFINE;
      }
    }
    else
    {
      cell_adapt_ops.assign(nb_cells, CellTransform::UNIFORM_REFINE);
    }

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

    mesh2d.adapt(adapt_schedule);
    dof_handler.adapt(adapt_schedule);
  }

  // --------------------------------------------------------------------------

  // Write the mesh in the output file
  const std::string outfilename2d = "h_adapted_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize(dof_handler.nb_active_cells());
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "active_cell_id");

  // --------------------------------------------------------------------------

  check_adapted_mesh(mesh2d, "geo_dofs", P4);

  const Real vol = tools::MeshInspector<Cart2D>::check_mesh_consistency(mesh2d, "geo_dofs", P4);
  std::cout << "Total evaluated domain volume = " << vol << std::endl;
  BOOST_CHECK_CLOSE(vol, 1.0, 1.e-10);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(circle_refinement_utest)
{
  std::cout << "\n************ RUNNING CIRCLE MESH h-ADAPTATION TEST - 2D "
               "************"
            << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));

  MeshAdaptSequence<Cart2D> adapt_schedule;

  // ---------
  std::vector<Uint> colors;
  interpolation::ScalarMeshFunction<Uint> marked_cells("", "marked_cells");
  interpolation::ScalarMeshFunction<Uint> cell_colors("", "cell_colors");
  // ---------

  const Uint nb_adapt_passes = 6;
  std::vector<CellTransform> cell_adapt_ops;

  const std::clock_t cpu_start = std::clock();

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    std::cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << std::endl;
    std::cout << "Adaptation pass #" << p + 1 << std::endl;
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    for (Uint c = 0; c < nb_cells; ++c)
    {
      const CellTopologyView<Cart2D> cell_tview    = mesh2d.active_cell(ActiveIdx(c));
      const CellGeometry<Cart2D::GDIM> cell_coords = cell_tview.coordinates();
      if (circle_intersects_element(0.5, 0.5, 0.3, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
      }
    }

    // ----------------------------------------------------
    // Write the coloring to a file (for debugging)

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

    const bool write_colors_to_file = false;
    if (write_colors_to_file)
    {
      colors.resize(dof_handler.nb_active_cells());

      mesh2d.color_cells_red_green(cell_adapt_ops, colors);

      marked_cells.resize(dof_handler.nb_active_cells());

      for (Uint i = 0; i < cell_adapt_ops.size(); ++i)
      {
        if (cell_adapt_ops[i] != CellTransform::NO_TRANS)
        {
          marked_cells[i] = 1;
        }
        else
        {
          marked_cells[i] = 0;
        }
      }

      cell_colors.resize(colors.size());
      for (Uint i = 0; i < colors.size(); ++i)
      {
        cell_colors[i] = colors[i];
      }

      gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "red_green.msh");
      gmshwriter.append_cell_function_to_file(mesh2d, "red_green.msh", marked_cells,
                                              "adaptation marker", 0, 0);
      gmshwriter.append_cell_function_to_file(mesh2d, "red_green.msh", cell_colors, "cell_colors",
                                              0, 0);
    }

    mesh2d.adapt(adapt_schedule);
    // dof_handler.adapt_update(mesh2d);
    dof_handler.adapt(adapt_schedule);
  }

  const std::clock_t cpu_end = std::clock();
  const double cpu_duration  = (cpu_end - cpu_start) / (double)CLOCKS_PER_SEC;

  std::cout.precision(10);
  std::cout << "\nMesh adaptation (" << nb_adapt_passes << " passes) took " << cpu_duration << "s\n"
            << std::endl;

  // Write the mesh in the output file
  const std::string outfilename2d = "circle_adapted_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize(dof_handler.nb_active_cells());
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label, "active_cell_id");

  // Finally, save the refinement level of each cell
  for (Uint i = 0; i < cell_label.nb_entries(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(mesh::ActiveIdx(i));
    cell_label[i]                      = topo_cell.refinement_level();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label,
                                          "active_cell_refinement_level");

  // --------------------------------------------------------------------------

  check_adapted_mesh(mesh2d, "geo_dofs", P4);

  const Real vol = tools::MeshInspector<Cart2D>::check_mesh_consistency(mesh2d, "geo_dofs", P4);
  std::cout << "Total evaluated domain volume (circle adapted mesh) = " << vol << std::endl;
  BOOST_CHECK_CLOSE(vol, 1.0, 1.e-10);

  mesh2d.write_dual_graph_to_gmsh("geo_dofs", "circle_adapted_dual_graph.msh");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(circle_mesh_coarsening_utest)
{
  std::cout << "\n************ RUNNING CIRCLE MESH COARSENING TEST - 2D "
               "************"
            << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));

  MeshAdaptSequence<Cart2D> adapt_schedule;

  const Uint nb_adapt_passes = 6;
  std::vector<CellTransform> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    for (Uint c = 0; c < nb_cells; ++c)
    {
      const CellTopologyView<Cart2D> cell_tview    = mesh2d.active_cell(ActiveIdx(c));
      const CellGeometry<Cart2D::GDIM> cell_coords = cell_tview.coordinates();
      if (circle_intersects_element(0.5, 0.5, 0.3, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
      }
    }

    // ----------------------------------------------------
    // Write the coloring to a file (for debugging)

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

    mesh2d.adapt(adapt_schedule);
    // dof_handler.adapt_update(mesh2d);
    dof_handler.adapt(adapt_schedule);
  }

  for (Uint cycle = 0; cycle < nb_adapt_passes - 2; ++cycle)
  // for (Uint cycle = 0; cycle < 1; ++cycle)
  {
    cell_adapt_ops.resize(mesh2d.nb_active_cells());
    cell_adapt_ops.assign(mesh2d.nb_active_cells(), CellTransform::NO_TRANS);

    for (Uint ac = 0; ac < mesh2d.nb_active_cells(); ++ac)
    {
      const mesh::CellTopologyView<Cart2D> tcell = mesh2d.active_cell(mesh::ActiveIdx(ac));
      if (tcell.refinement_level() == nb_adapt_passes - cycle)
      {
        cell_adapt_ops[ac] = CellTransform::COARSEN;
      }
    }

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::coarsen);

    mesh2d.adapt(adapt_schedule);
    dof_handler.adapt_update(mesh2d);
  }
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "circle_coarsened.msh");
  mesh2d.write_dual_graph_to_gmsh("geo_dofs", "circle_coarsened_dual.msh");
}

// ----------------------------------------------------------------------------

// The purpose of this test is to verify that node-to-cell
// inicidence computation works correctly on adapted meshes.

BOOST_AUTO_TEST_CASE(node_to_cell_comp_in_adapted_mesh)
{
  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri_no_bdry.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  result_of::dof_map_t<Cart2D> &dof_handler = *(mesh2d.dof_storage("geo_dofs"));

  MeshAdaptSequence<Cart2D> adapt_schedule;

  const Uint nb_adapt_passes = 8;
  std::vector<CellTransform> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    for (Uint c = 0; c < nb_cells; ++c)
    {
      const CellTopologyView<Cart2D> cell_tview    = mesh2d.active_cell(ActiveIdx(c));
      const CellGeometry<Cart2D::GDIM> cell_coords = cell_tview.coordinates();
      if (circle_intersects_element(0.5, 0.5, 0.3, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
      }
    }

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

    mesh2d.adapt(adapt_schedule);
    dof_handler.adapt(adapt_schedule);
  }

  // Write the mesh in the output file
  const std::string outfilename2d = "node_to_cell_adapted_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  // Compute node-to-cell incidences
  common::BlockArray<std::tuple<Uint, Uint>, Uint> incident_p1_nodes;
  common::BlockArray<Uint, Uint> cell_to_p1_node;
  mesh2d.compute_node_to_cell_connectivity(incident_p1_nodes, cell_to_p1_node);

  /*
  for (Uint b = 0; b < incident_p1_nodes.nb_blocks(); ++b)
  {
    const common::ArrayView<const std::tuple<Uint, Uint>, Uint> block =
  incident_p1_nodes.const_block(b);
    for (Uint i = 0; i < block.size(); ++i)
    {
      std::cout << "[" << std::get<0>(block[i]) + 1 << "," <<
  std::get<1>(block[i]) << "] ";
    }
    std::cout << std::endl;
  }
  */

  adapt::LocalInterpolator loc_interpolator;
  std::vector<Real> ref_node_coord;

  for (Uint b = 0; b < incident_p1_nodes.nb_blocks(); ++b)
  {
    common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> incident_nodes =
        incident_p1_nodes.const_block(b);

    BOOST_CHECK_GE(incident_nodes.size(), 1);

    /*
    const Uint active_cell_id = std::get<0>(incident_nodes[0]);
    const Uint local_node_id = std::get<1>(incident_nodes[0]);

    const mesh::MeshEntity active_cell =
    dof_handler.active_cell(active_cell_id); std::cout << active_cell << "
    "; std::cout << "[" << active_cell_id + 1 << " -> " <<
    active_cell.vertex(local_node_id) + 1
              << "]";
    */

    const ActiveIdx ref_cell_active_idx(std::get<0>(incident_nodes[0]));
    const CellTopologyView<Cart2D> ref_cell_tview = dof_handler.tcell(ref_cell_active_idx);
    const MeshEntity ref_active_cell              = dof_handler.active_cell(ref_cell_active_idx);

    const math::DenseConstMatView<Real> ref_cell_coords = loc_interpolator.transfer_coords(
        ref_cell_tview.pt_set_id(), ref_active_cell.pt_set_id(), ref_cell_tview.coordinates());

    const math::DenseConstVecView<Real> tmp =
        ref_cell_coords.row_transpose(std::get<1>(incident_nodes[0]));

    ref_node_coord.resize(tmp.size());
    for (Uint i = 0; i < tmp.size(); ++i)
    {
      ref_node_coord[i] = tmp[i];
    }

    for (Uint n = 1; n < incident_nodes.size(); ++n)
    {
      /*
      const Uint active_cell_id = std::get<0>(incident_nodes[n]);
      const Uint local_node_id = std::get<1>(incident_nodes[n]);
      const mesh::MeshEntity active_cell =
      dof_handler.active_cell(active_cell_id);

      std::cout << "    " << active_cell << " ";
      std::cout << " [" << active_cell_id + 1 << " -> " <<
      active_cell.vertex(local_node_id) + 1
                << "]";
      */

      const ActiveIdx cell_active_idx(std::get<0>(incident_nodes[n]));
      const CellTopologyView<Cart2D> cell_tview = dof_handler.tcell(cell_active_idx);
      const MeshEntity active_cell              = dof_handler.active_cell(cell_active_idx);

      const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
          cell_tview.pt_set_id(), active_cell.pt_set_id(), cell_tview.coordinates());

      const math::DenseConstVecView<Real> node_coord =
          ref_cell_coords.row_transpose(std::get<1>(incident_nodes[n]));

      const Real dx0  = node_coord[X0] - ref_node_coord[X0];
      const Real dx1  = node_coord[X1] - ref_node_coord[X1];
      const Real dist = std::sqrt(dx0 * dx0 + dx1 * dx1);

      BOOST_CHECK_LE(dist, 1.e-10);
    }
    // std::cout << std::endl;
  }

  // Check that the arrays incident_p1_nodes and cell_to_p1_node are
  // consistent
  for (Uint cell_id = 0; cell_id < cell_to_p1_node.nb_blocks(); ++cell_id)
  {
    // std::cout << "Active cell [" << cell_id << "] has p1 vertices in
    // blocks:"
    // << std::endl;
    const common::ArrayView<const Uint, _1D, Uint> p1_vert_pos =
        cell_to_p1_node.const_block(cell_id);

    for (Uint local_p1_node_idx = 0; local_p1_node_idx < p1_vert_pos.size(); ++local_p1_node_idx)
    {
      bool node_found_in_block = false;
      const common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> incident_vert_block =
          incident_p1_nodes.const_block(p1_vert_pos[local_p1_node_idx]);
      // std::cout << "   { ";
      for (Uint j = 0; j < incident_vert_block.size(); ++j)
      {
        // std::cout << "(" << std::get<0>(incident_vert_block[j]) <<
        // ","
        //           << std::get<1>(incident_vert_block[j]) << ") ";
        if ((std::get<0>(incident_vert_block[j]) == cell_id) &&
            (std::get<1>(incident_vert_block[j]) == local_p1_node_idx))
        {
          node_found_in_block = true;
        }
      }
      BOOST_CHECK_EQUAL(node_found_in_block, true);
      // std::cout << " }" << std::endl;
    }
    // std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
