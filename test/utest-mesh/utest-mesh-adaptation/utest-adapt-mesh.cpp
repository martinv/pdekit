/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE adapt_mesh_test
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>

#include "common/Constants.hpp"
#include "common/StringUtils.hpp"
#include "mesh/DofCoordinates.hpp"
#include "mesh/adaptation/GeometryAdapter.hpp"
#include "mesh/adaptation/MeshAdaptSequence.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct AdaptMeshUtestFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(AdaptMesh_TestSuite, AdaptMeshUtestFixture)

// ----------------------------------------------------------------------------
#if 1
BOOST_AUTO_TEST_CASE(adapt_cell_geometry)
{
  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri_no_bdry.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");

  const CellTopologyView<Cart2D> tview   = (*dof_handler).tcell(ActiveIdx(27));
  CellGeometry<Cart2D::GDIM> cell_coords = tview.coordinates();

  adapt::GeometryAdapter<Cart2D> geo_adapter;

  adapt::CellAdaptOpTag split_tag(ElemShape::Triag, CellTransform::UNIFORM_REFINE);

  const common::ArrayView<const math::DenseDMat<Real>, _1D, Uint> child_coords =
      geo_adapter.compute_child_coords(tview.pt_set_id(), split_tag, cell_coords);

  /*
  std::vector<Real> cell_coord_vec(cell_coords.size() * Cart2D::GDIM);
  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    const math::DenseConstVecView<Real> one_node_coords = cell_coords.c(n);
    for (Uint d = 0; d < Cart2D::GDIM; ++d)
    {
      cell_coord_vec[n * Cart2D::GDIM + d] = one_node_coords[d];
    }
  }

  const math::DenseConstMatView<Real> cell_coord_mat(cell_coord_vec.data(),
  Cart2D::GDIM, cell_coord_vec.size(), Cart2D::GDIM);

  const common::ArrayView<const math::DenseDMat<Real> const, Uint>
  child_coords = geo_adapter.compute_child_coords(cell.std_region_id(),
  split_tag, cell_coord_mat);
  */

  // Write the coordinates into files - can be plotted with gnuplot
  std::ofstream outfile;
  outfile.precision(14);
  outfile.setf(std::ios::fixed);

  outfile.open("coords_parent_cell.dat");
  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    outfile << cell_coords.const_node_view(n) << std::endl;
  }
  outfile.close();

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    outfile.open("coords_child_cell" + common::StringUtils::to_string(c) + ".dat");
    outfile << child_coords[c];
    outfile.close();
  }

  // const std::string outfilename2d = "tmp_geo_adapt_" + infilename2d;
  // gmshwriter.write_mesh_to_file(mesh2d, outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(h_adapt_mesh_2D_interior_cell_hanging_nodes)
{
  std::cout << "\n************ RUNNING MESH h-ADAPTATION TEST (HANGING NODES, "
               "INTERIOR CELL) - 2D "
               "************"
            << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  /*
  for(Uint c = 0; c < dof_handler.nb_cells(); ++c)
  {
    std::cout << "[" << c << "] " << dof_handler.cell(c) << std::endl;
  }
  */

  Uint nb_cells = mesh2d.nb_active_cells();

  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;
  cell_adapt_ops[9] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2d.adapt(adapt_schedule);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  (*dof_handler).adapt(adapt_schedule);

  // Second adaptation pass
  /*
  nb_cells = mesh2d.nb_active_cells();
  cell_adapt_ops.resize(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::DO_NOTHING);
  cell_adapt_ops[17] = CellTransform::UNIFORM_REFINE;
  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops,
  h_AdaptStrategy::w_hanging_nodes); mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);
  */

  CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(5));

  std::cout << topo_cell << std::endl;

  // Write the mesh in the output file
  const std::string outfilename2d = "h_adapted_interior_cell_hanging_nodes_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  std::vector<Uint> cell_label;
  cell_label.resize((*dof_handler).nb_active_cells());
  const common::ArrayView<const Uint, _1D, Uint> cell_label_view(cell_label.data(),
                                                                 cell_label.size());

  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label_view,
                                          "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label_view, "active_cell_id");

  /*
  typedef result_of::dof_handler<Cart2D>::type::const_dof_iterator
  const_dof_iterator; for (const_dof_iterator dit = dof_handler.cbegin(); dit
  != dof_handler.cend(); ++dit)
  {
    const MeshEntity cell = *dit;
    std::cout << "[" << cell.idx() << "] " << cell << std::endl;
  }
  */

  mesh::CellPath path(FlatIdx(8), {0, 1, 1, 2, 1});
  const mesh::CellTopologyView<Cart2D> leaf_cell = mesh2d.path_leaf(path);
  std::cout << "Leaf cell = " << leaf_cell.active_idx() << std::endl;
}
#endif
// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(h_adapt_mesh_2D_interior_cell_red_green)
{
  std::cout << "\n************ RUNNING MESH h-ADAPTATION TEST (RED-GREEN, "
               "INTERIOR CELL) - 2D "
               "************"
            << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  /*
  for(Uint c = 0; c < dof_handler.nb_cells(); ++c)
  {
    std::cout << "[" << c << "] " << dof_handler.cell(c) << std::endl;
  }
  */

  Uint nb_cells = mesh2d.nb_active_cells();

  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;
  std::vector<CellTransform> cell_adapt_ops(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[8] = CellTransform::UNIFORM_REFINE;
  cell_adapt_ops[9] = CellTransform::UNIFORM_REFINE;
  // cell_adapt_ops[17] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);
  mesh2d.adapt(adapt_schedule);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  //(*dof_handler).adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);

  /*
  std::cout << "Printing skeleton:" << std::endl;
  mesh2d.print_complete_skeleton(_1D);
  */

  // Second adaptation pass
  nb_cells = mesh2d.nb_active_cells();
  cell_adapt_ops.resize(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[13] = CellTransform::UNIFORM_REFINE;
  // cell_adapt_ops[23] = CellTransform::UNIFORM_REFINE;
  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);
  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);

  // Third adaptation pass
  std::cout << "RED-GREEN PASS #3" << std::endl;

  nb_cells = mesh2d.nb_active_cells();
  cell_adapt_ops.resize(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[30] = CellTransform::UNIFORM_REFINE;
  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);

  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);

  // Fourth adaptation pass

  std::cout << "RED-GREEN PASS #4" << std::endl;

  nb_cells = mesh2d.nb_active_cells();
  cell_adapt_ops.resize(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::NO_TRANS);
  cell_adapt_ops[38] = CellTransform::UNIFORM_REFINE;
  cell_adapt_ops[41] = CellTransform::UNIFORM_REFINE;
  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::red_green);

  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);

  /*
  std::cout << "RED-GREEN PASS #4" << std::endl;
  nb_cells = mesh2d.nb_active_cells();
  cell_adapt_ops.resize(nb_cells);
  cell_adapt_ops.assign(nb_cells, CellTransform::DO_NOTHING);


  cell_adapt_ops[24] = CellTransform::COARSEN;
  cell_adapt_ops[25] = CellTransform::COARSEN;
  cell_adapt_ops[26] = CellTransform::COARSEN;
  cell_adapt_ops[27] = CellTransform::COARSEN;
  cell_adapt_ops[28] = CellTransform::COARSEN;
  cell_adapt_ops[29] = CellTransform::COARSEN;
  cell_adapt_ops[30] = CellTransform::COARSEN;
  cell_adapt_ops[31] = CellTransform::COARSEN;
  cell_adapt_ops[32] = CellTransform::COARSEN;
  cell_adapt_ops[33] = CellTransform::COARSEN;
  cell_adapt_ops[34] = CellTransform::COARSEN;
  cell_adapt_ops[35] = CellTransform::COARSEN;
  cell_adapt_ops[36] = CellTransform::COARSEN;
  cell_adapt_ops[37] = CellTransform::COARSEN;
  cell_adapt_ops[38] = CellTransform::COARSEN;
  cell_adapt_ops[39] = CellTransform::COARSEN;
  */

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::coarsen);

  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);

  CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(5));

  std::cout << topo_cell << std::endl;

  // Write the mesh in the output file
  const std::string outfilename2d = "h_adapted_interior_cell_red_green_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  std::vector<Uint> cell_label;
  cell_label.resize((*dof_handler).nb_active_cells());
  const common::ArrayView<const Uint, _1D, Uint> cell_label_view(cell_label.data(),
                                                                 cell_label.size());

  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label_view,
                                          "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, outfilename2d, cell_label_view, "active_cell_id");
  mesh2d.write_dual_graph_to_gmsh("geo_dofs", "dual_" + outfilename2d);
}

// ----------------------------------------------------------------------------

#if 1
BOOST_AUTO_TEST_CASE(coarsen_mesh_2D_interior_cell)
{
  std::cout << "\n************ RUNNING MESH COARSENING TEST (INTERIOR CELL) - "
               "2D ************"
            << std::endl;

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
  cell_adapt_ops[9] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2d.adapt(adapt_schedule);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  (*dof_handler).adapt(adapt_schedule);

  // Write the mesh in the output file
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "coarsened_interior_cell_before.msh");

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  std::vector<Uint> cell_label;
  cell_label.resize((*dof_handler).nb_active_cells());
  const common::ArrayView<const Uint, _1D, Uint> cell_label_view(cell_label.data(),
                                                                 cell_label.size());

  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, "coarsened_interior_cell_before.msh",
                                          cell_label_view, "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, "coarsened_interior_cell_before.msh",
                                          cell_label_view, "active_cell_id");

  std::cout << "************************************************" << std::endl;

  cell_adapt_ops.resize(mesh2d.nb_active_cells());
  cell_adapt_ops.assign(mesh2d.nb_active_cells(), CellTransform::NO_TRANS);
  cell_adapt_ops[19] = CellTransform::COARSEN;
  // cell_adapt_ops[22] = CellTransform::COARSEN;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::coarsen);

  adapt_schedule.print_ops();

  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "coarsened_interior_cell_after.msh");
  mesh2d.write_dual_graph_to_gmsh("geo_dofs", "coarsened_interior_cell_after_dual.msh");
  std::cout << "************************************************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(h_adapt_mesh_2D_bdry_cell)
{
  std::cout << "\n************ RUNNING MESH h-ADAPTATION TEST (BOUNDARY CELL) "
               "- 2D ************"
            << std::endl;

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
  cell_adapt_ops[13] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2d.adapt(adapt_schedule);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  (*dof_handler).adapt(adapt_schedule);

  CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(5));

  std::cout << topo_cell << std::endl;

  // Write the mesh in the output file
  const std::string outfilename2d = "h_adapted_bdry_cell_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(coarsen_mesh_2D_bdry_cell)
{
  std::cout << "\n************ RUNNING MESH COARSENING TEST (BOUNDARY CELL) - "
               "2D ************"
            << std::endl;

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
  cell_adapt_ops[13] = CellTransform::UNIFORM_REFINE;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::w_hanging_nodes);
  mesh2d.adapt(adapt_schedule);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  (*dof_handler).adapt(adapt_schedule);

  // Write the mesh in the output file
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "coarsened_boundary_cell_before.msh");

  // Make a function that displays active cell numbers as they are numbered
  // inside the code
  std::vector<Uint> cell_label;
  cell_label.resize((*dof_handler).nb_active_cells());
  const common::ArrayView<const Uint, _1D, Uint> cell_label_view(cell_label.data(),
                                                                 cell_label.size());

  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.linear_pos_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, "coarsened_boundary_cell_before.msh",
                                          cell_label_view, "absolute_cell_id");

  // Now reuse the same function to save active cell ids
  for (Uint i = 0; i < cell_label.size(); ++i)
  {
    CellTopologyView<Cart2D> topo_cell = mesh2d.active_cell(ActiveIdx(i));
    cell_label[i]                      = topo_cell.active_idx().id();
  }
  gmshwriter.append_cell_function_to_file(mesh2d, "coarsened_boundary_cell_before.msh",
                                          cell_label_view, "active_cell_id");

  std::cout << "************************************************" << std::endl;

  cell_adapt_ops.resize(mesh2d.nb_active_cells());
  cell_adapt_ops.assign(mesh2d.nb_active_cells(), CellTransform::NO_TRANS);
  cell_adapt_ops[20] = CellTransform::COARSEN;
  // cell_adapt_ops[22] = CellTransform::COARSEN;

  adapt_schedule.define_h_adapt_ops(mesh2d, cell_adapt_ops, h_AdaptStrategy::coarsen);

  adapt_schedule.print_ops();

  mesh2d.adapt(adapt_schedule);
  (*dof_handler).adapt_update(mesh2d);
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", "coarsened_boundary_cell_after.msh");
  std::cout << "************************************************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(p_adapt_mesh_2D)
{
  std::cout << "\n************ RUNNING MESH p-ADAPTATION TEST - 2D ************" << std::endl;

  // Create a mesh
  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "unit_square_mini.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  /*
  for(Uint c = 0; c < dof_handler.nb_cells(); ++c)
  {
    std::cout << "[" << c << "] " << dof_handler.cell(c) << std::endl;
  }
  */

  const Uint nb_cells = mesh2d.nb_active_cells();
  std::vector<Uint> cell_p_order(nb_cells, P1);
  cell_p_order[8] = P2;
  cell_p_order[9] = P3;

  // First adaptation pass
  MeshAdaptSequence<Cart2D> adapt_schedule;
  adapt_schedule.define_p_adapt_ops(cell_p_order);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> dof_handler = mesh2d.dof_storage("geo_dofs");
  (*dof_handler).adapt(adapt_schedule);

  // Second adaptation pass
  cell_p_order[8] = P5;
  adapt_schedule.define_p_adapt_ops(cell_p_order);
  (*dof_handler).adapt(adapt_schedule);

  // Write the mesh in the output file
  const std::string outfilename2d = "p_adapted_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
}
#endif
// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
