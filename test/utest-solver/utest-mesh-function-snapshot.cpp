/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_function_snapshot_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "interpolation/FunctionSpace.hpp"
#include "interpolation/mesh_function/MeshFunctionSnapshot.hpp"
#include "interpolation/mesh_function/ScalarMeshFunction.hpp"
#include "interpolation/mesh_function/VectorMeshFunction.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

template <Uint GeoDim>
bool circle_intersects_element(const Real xc, const Real yc, const Real R,
                               const CellGeometry<GeoDim> &cell_coords)
{
  Uint nb_nodes_in  = 0;
  Uint nb_nodes_out = 0;

  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    const math::DenseVecView<const Real> point_coords = cell_coords.const_node_view(n);
    const Real node_dist = std::sqrt((point_coords[X0] - xc) * (point_coords[X0] - xc) +
                                     (point_coords[X1] - yc) * (point_coords[X1] - yc));

    if (node_dist < R)
    {
      nb_nodes_in++;
    }
    else
    {
      nb_nodes_out++;
    }
    if (nb_nodes_in * nb_nodes_out > 0)
    {
      return true;
    }
  }
  return (nb_nodes_in * nb_nodes_out > 0);
}

// ----------------------------------------------------------------------------

struct MeshFunctSnapshotUtestFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MeshFunctionSnapshot_TestSuite, MeshFunctSnapshotUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_mesh_function_snapshot_h_adapt_utest)
{
  typedef mesh::Tria<Cart2D> mesh_type;
  typedef result_of::dof_map_t<Cart2D> dof_map_type;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  MeshCreator::make_unit_quad(mesh2d, "geo_dofs", 20, true);

  common::PtrHandle<dof_map_type> geo_dof_handler = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dof_handler = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_discontinuous(mesh2d, *geo_dof_handler, *sol_dof_handler, P3,
                                    PointSetID::Warpblend);

  interpolation::ScalarMeshFunction<Real> u_scalar("", "u_scalar");
  u_scalar.resize((*sol_dof_handler).nb_nodes());

  interpolation::ScalarMeshFunction<Uint> u_scalar_cellwise("", "u_scalar_cellwise");
  u_scalar_cellwise.resize((*sol_dof_handler).nb_active_cells());

  interpolation::VectorMeshFunction<Real> u_vector("", "u_vector");
  u_vector.resize(3, (*sol_dof_handler).nb_nodes());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < (*sol_dof_handler).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<Cart2D> tcell_view = (*sol_dof_handler).tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell = (*sol_dof_handler).active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const math::DenseConstVecView<Real> node_coords = cell_coords.row_transpose(v);
      u_scalar[cell.vertex(v)] = std::sin(4. * node_coords[X0]) * std::cos(3. * node_coords[X1]);

      interpolation::VectorMeshFunction<Real>::entry_type node_value =
          u_vector.value(cell.vertex(v));

      node_value[0] = std::sin(2. * node_coords[X0]) * std::cos(3. * node_coords[X1]);
      node_value[1] = std::cos(4. * node_coords[X0]) * std::sin(3. * node_coords[X1]);
      node_value[2] = std::cos(5. * node_coords[X0]) * std::cos(6. * node_coords[X1]);
    }
  }

  MeshAdaptSequence<Cart2D> adapt_schedule;
  interpolation::MeshFunctionSnapshot<Real> u_snapshot;

  const Uint nb_adapt_passes             = 6;
  const Uint cellwise_snapshot_threshold = 3;
  std::vector<CellTransform> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    // u_snapshot.clear();
    u_snapshot.create<Cart2D>(mesh2d, *sol_dof_handler, u_scalar);
    u_snapshot.create<Cart2D>(mesh2d, *sol_dof_handler, u_vector);

    std::cout << "Adaptation pass #" << p + 1 << std::endl;
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    if (p == cellwise_snapshot_threshold)
    {
      u_scalar_cellwise.resize(mesh2d.nb_active_cells());
      u_scalar_cellwise.fill(0);
    }

    for (Uint c = 0; c < nb_cells; ++c)
    {
      const mesh::CellTopologyView<Cart2D> tcell_view =
          (*sol_dof_handler).tcell(mesh::ActiveIdx(c));
      const mesh::CellGeometry<Cart2D::GDIM> cell_coords = tcell_view.coordinates();

      if (circle_intersects_element(0.5, 0.5, 0.3, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
        if (p == cellwise_snapshot_threshold)
        {
          u_scalar_cellwise[c] = 1;
        }
      }

      if (circle_intersects_element(-0.5, -0.5, 0.25, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
        if (p == cellwise_snapshot_threshold)
        {
          u_scalar_cellwise[c] = 1;
        }
      }
    }

    if (p >= cellwise_snapshot_threshold)
    {
      u_snapshot.create_cellwise<Cart2D>(mesh2d, *sol_dof_handler, u_scalar_cellwise);
    }

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);

    mesh2d.adapt(adapt_schedule);
    (*geo_dof_handler).adapt(adapt_schedule);
    (*sol_dof_handler).adapt(adapt_schedule);

    u_snapshot.restore_function<Cart2D>(mesh2d, *sol_dof_handler, u_scalar);
    u_snapshot.restore_function<Cart2D>(mesh2d, *sol_dof_handler, u_vector);

    if (p >= cellwise_snapshot_threshold)
    {
      u_snapshot.restore_function_cellwise<Cart2D>(mesh2d, *sol_dof_handler, u_scalar_cellwise);
    }
  }

  // Write the mesh in the output file
  const std::string outfilename2d = "mesh_function_h_adapted.msh";
  gmshwriter.write_mesh_to_file(mesh2d, "sol_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize((*geo_dof_handler).nb_active_cells());
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
  gmshwriter.append_nodal_function_to_file(mesh2d, outfilename2d, u_scalar, "u_scalar");
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, u_scalar_cellwise,
                                          "u_scalar_cellwise");
  gmshwriter.append_nodal_function_to_file(mesh2d, outfilename2d, u_vector, "u_vector");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cellwise_scalar_mesh_function_snapshot_h_adapt_utest)
{
  typedef mesh::Tria<Cart2D> mesh_type;
  typedef result_of::dof_map_t<Cart2D> dof_map_type;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  MeshCreator::make_unit_quad(mesh2d, "geo_dofs", 20, true);

  common::PtrHandle<dof_map_type> geo_dof_handler = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dof_handler = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_discontinuous(mesh2d, *geo_dof_handler, *sol_dof_handler, P3,
                                    PointSetID::Warpblend);

  interpolation::ScalarMeshFunction<Uint> u_scalar_cellwise("", "u_scalar_cellwise");
  u_scalar_cellwise.resize((*sol_dof_handler).nb_active_cells());

  MeshAdaptSequence<Cart2D> adapt_schedule;
  interpolation::MeshFunctionSnapshot<Uint> u_snapshot;

  const Uint nb_adapt_passes             = 6;
  const Uint cellwise_snapshot_threshold = 3;
  std::vector<CellTransform> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    // u_snapshot.clear();

    std::cout << "Adaptation pass #" << p + 1 << std::endl;
    const Uint nb_cells = mesh2d.nb_active_cells();

    // cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
    if (p == cellwise_snapshot_threshold)
    {
      u_scalar_cellwise.resize(mesh2d.nb_active_cells());
      u_scalar_cellwise.fill(0);
    }

    for (Uint c = 0; c < nb_cells; ++c)
    {
      const mesh::CellTopologyView<Cart2D> tcell_view =
          (*sol_dof_handler).tcell(mesh::ActiveIdx(c));
      const mesh::CellGeometry<Cart2D::GDIM> cell_coords = tcell_view.coordinates();

      if (circle_intersects_element(0.5, 0.5, 0.3, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
        if (p == cellwise_snapshot_threshold)
        {
          u_scalar_cellwise[c] = 1;
        }
      }

      if (circle_intersects_element(-0.5, -0.5, 0.25, cell_coords))
      {
        cell_adapt_ops[c] = CellTransform::UNIFORM_REFINE;
        if (p == cellwise_snapshot_threshold)
        {
          u_scalar_cellwise[c] = 1;
        }
      }
    }

    if (p >= cellwise_snapshot_threshold)
    {
      u_snapshot.create_cellwise<Cart2D>(mesh2d, *sol_dof_handler, u_scalar_cellwise);
    }

    adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);

    mesh2d.adapt(adapt_schedule);
    (*geo_dof_handler).adapt_update(mesh2d);
    (*sol_dof_handler).adapt_update(mesh2d);

    if (p >= cellwise_snapshot_threshold)
    {
      u_snapshot.restore_function_cellwise<Cart2D>(mesh2d, *sol_dof_handler, u_scalar_cellwise);
    }
  }

  // Write the mesh in the output file
  const std::string outfilename2d = "cellwise_mesh_function_h_adapted.msh";
  gmshwriter.write_mesh_to_file(mesh2d, "sol_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize((*geo_dof_handler).nb_active_cells());
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

  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, u_scalar_cellwise,
                                          "u_scalar_cellwise");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scalar_mesh_function_snapshot_p_adapt_utest)
{
  typedef mesh::Tria<Cart2D> mesh_type;
  typedef result_of::dof_map_t<Cart2D> dof_map_type;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  MeshCreator::make_unit_quad(mesh2d, "geo_dofs", 20, true);

  common::PtrHandle<dof_map_type> geo_dof_handler = mesh2d.dof_storage("geo_dofs");
  common::PtrHandle<dof_map_type> sol_dof_handler = mesh2d.create_dof_storage("sol_dofs");

  dof_map_type::clone_discontinuous(mesh2d, *geo_dof_handler, *sol_dof_handler, P2,
                                    PointSetID::Warpblend);

  interpolation::ScalarMeshFunction<Real> u_scalar("", "u_scalar");
  u_scalar.resize((*sol_dof_handler).nb_nodes());

  interpolation::VectorMeshFunction<Real> u_vector("", "u_vector");
  u_vector.resize(3, (*sol_dof_handler).nb_nodes());

  mesh::adapt::LocalInterpolator loc_interpolator;

  for (Uint c = 0; c < (*sol_dof_handler).nb_active_cells(); ++c)
  {
    const mesh::CellTopologyView<Cart2D> tcell_view = (*sol_dof_handler).tcell(mesh::ActiveIdx(c));
    const mesh::MeshEntity cell = (*sol_dof_handler).active_cell(mesh::ActiveIdx(c));

    const math::DenseConstMatView<Real> cell_coords = loc_interpolator.transfer_coords(
        tcell_view.pt_set_id(), cell.pt_set_id(), tcell_view.coordinates());

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const math::DenseConstVecView<Real> node_coords = cell_coords.row_transpose(v);
      u_scalar[cell.vertex(v)] = std::sin(4. * node_coords[X0]) * std::cos(3. * node_coords[X1]);

      interpolation::VectorMeshFunction<Real>::entry_type node_value =
          u_vector.value(cell.vertex(v));

      node_value[0] = std::sin(2. * node_coords[X0]) * std::cos(3. * node_coords[X1]);
      node_value[1] = std::cos(4. * node_coords[X0]) * std::sin(3. * node_coords[X1]);
      node_value[2] = std::cos(5. * node_coords[X0]) * std::cos(6. * node_coords[X1]);
    }
  }

  MeshAdaptSequence<Cart2D> adapt_schedule;
  interpolation::MeshFunctionSnapshot<Real> u_snapshot;

  const Uint nb_adapt_passes = 3;
  std::vector<Uint> cell_adapt_ops;

  for (Uint p = 0; p < nb_adapt_passes; ++p)
  {
    // u_snapshot.clear();
    u_snapshot.create<Cart2D>(mesh2d, *sol_dof_handler, u_scalar);
    u_snapshot.create<Cart2D>(mesh2d, *sol_dof_handler, u_vector);

    std::cout << "Adaptation pass #" << p + 1 << std::endl;
    const Uint nb_cells = mesh2d.nb_active_cells();

    cell_adapt_ops.resize(nb_cells);
    cell_adapt_ops.assign(nb_cells, 3 + p);

    adapt_schedule.define_p_adapt_ops(cell_adapt_ops);

    // mesh2d.adapt(adapt_schedule);
    (*geo_dof_handler).adapt(adapt_schedule);
    (*sol_dof_handler).adapt(adapt_schedule);

    u_snapshot.restore_function<Cart2D>(mesh2d, *sol_dof_handler, u_scalar);
    u_snapshot.restore_function<Cart2D>(mesh2d, *sol_dof_handler, u_vector);
  }

  // Write the mesh in the output file
  const std::string outfilename2d = "mesh_function_p_adapted.msh";
  gmshwriter.write_mesh_to_file(mesh2d, "sol_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  interpolation::ScalarMeshFunction<Uint> cell_label("", "cell_id");
  cell_label.resize((*geo_dof_handler).nb_active_cells());
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
  gmshwriter.append_nodal_function_to_file(mesh2d, outfilename2d, u_scalar, "u_scalar");
  gmshwriter.append_nodal_function_to_file(mesh2d, outfilename2d, u_vector, "u_vector");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(local_elem_p_interpolation_utest)
{
  interpolation::detail::LocalElemPInterpolation<Real> local_elem_p_interp;

  const mesh::PointSetTag tag_elem_in(ElemShape::Quad, PolyOrderID::P2, PointSetID::Equidist);
  const mesh::PointSetTag tag_elem_out(ElemShape::Quad, PolyOrderID::P3, PointSetID::Equidist);

  math::DenseDMat<Real> values_in;
  mesh::StdRegion elem_in;

  elem_in.change_type(tag_elem_in);

  math::DenseDMat<Real> const &ref_coords = elem_in.get().coordinates();

  values_in.resize(elem_in.get().nb_nodes(), 3);

  for (Uint n = 0; n < elem_in.get().nb_nodes(); ++n)
  {
    const math::DenseConstVecView<Real> ref_coord = ref_coords.const_row_transp(n);
    values_in(n, 0)                               = ref_coord[X0] + ref_coord[X1];
    values_in(n, 1)                               = ref_coord[X0] * ref_coord[X0] + ref_coord[X1];
    values_in(n, 2)                               = ref_coord[X0] + ref_coord[X1] * ref_coord[X1];
  }

  local_elem_p_interp.interpolate(tag_elem_in, tag_elem_out, values_in);

  std::cout << "Values in = " << std::endl << values_in << std::endl;
  std::cout << "Interpolated values = " << std::endl
            << local_elem_p_interp.interpolated_values() << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
