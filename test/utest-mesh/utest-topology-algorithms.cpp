/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE topology_algorithms_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/Constants.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

struct TopologyAlgorithmsUtestFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(TopologyAlgorithms_TestSuite, TopologyAlgorithmsUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(discontinuous_topology_creation_2D_adaptive)
{
  std::cout << std::endl
            << "************ RUNNING DISCONTINUOUS CONNECTIVITY "
            << "CREATION TEST - 2D ADAPTIVE ************" << std::endl;

  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri.msh";

  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "continuous_dofs");

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> continuous_dofs =
      mesh2d.dof_storage("continuous_dofs");
  common::PtrHandle<Tria<Cart2D>::dof_storage_type> discontinuous_dofs =
      mesh2d.create_dof_storage("discontinuous_dofs");

  // Perform the actual upgrade
  Tria<Cart2D>::dof_storage_type::clone_discontinuous(mesh2d, *continuous_dofs, *discontinuous_dofs,
                                                      P5, PointSetID::Warpblend);

  // Write the mesh in the output file
  const std::string outfilename2d = "output_discontinuous_adaptive_" + infilename2d;

  gmshwriter.write_mesh_to_file(mesh2d, "discontinuous_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(discontinuous_topology_creation_3D_adaptive)
{
  std::cout << std::endl
            << "************ RUNNING DISCONTINUOUS CONNECTIVITY "
            << "CREATION TEST - 3D ADAPTIVE ************" << std::endl;

  Tria<Cart3D> mesh3d("mesh3d");

  // Read some data in it
  const std::string infilename3d = "cube.msh";

  gmshreader.read_mesh_from_file(infilename3d, mesh3d, "continuous_dofs");

  // Perform the actual upgrade
  common::PtrHandle<Tria<Cart3D>::dof_storage_type> continuous_dofs =
      mesh3d.dof_storage("continuous_dofs");
  common::PtrHandle<Tria<Cart3D>::dof_storage_type> discontinuous_dofs =
      mesh3d.create_dof_storage("discontinuous_dofs");

  // Perform the actual upgrade
  Tria<Cart3D>::dof_storage_type::clone_discontinuous(mesh3d, *continuous_dofs, *discontinuous_dofs,
                                                      P3, PointSetID::Equidist);

  // Write the mesh in the output file
  const std::string outfilename3d = "output_discontinuous_adaptive_" + infilename3d;

  gmshwriter.write_mesh_to_file(mesh3d, "discontinuous_dofs", outfilename3d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(node_to_element_connectivity_2D)
{
  Tria<Cart2D> mesh2d("mesh2d");
  common::BlockArray<std::tuple<Uint, Uint>, Uint> node_to_cell;
  common::BlockArray<Uint, Uint> cell_to_p1_node;

  // Read some data in it
  const std::string infilename2d = "test_p1_tri_no_bdry.msh";

  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "continuous_dofs");
  mesh2d.compute_node_to_cell_connectivity(node_to_cell, cell_to_p1_node);

  common::PtrHandle<Tria<Cart2D>::dof_storage_type> geo_dofs =
      mesh2d.dof_storage("continuous_dofs");
  for (Uint b = 0; b < node_to_cell.nb_blocks(); ++b)
  {
    common::ArrayView<const std::tuple<Uint, Uint>, _1D, Uint> incident_nodes =
        node_to_cell.const_block(b);

    BOOST_CHECK_GE(incident_nodes.size(), 1);

    const mesh::CellTopologyView<Cart2D> tcell_view_ref =
        (*geo_dofs).tcell(ActiveIdx(std::get<0>(incident_nodes[0])));
    const mesh::CellGeometry<Cart2D::GDIM> cell_coords_ref = tcell_view_ref.coordinates();

    const auto node_coord_ref = cell_coords_ref.const_node_view(std::get<1>(incident_nodes[0]));

    for (Uint n = 1; n < incident_nodes.size(); ++n)
    {
      const mesh::CellTopologyView<Cart2D> tcell_view =
          (*geo_dofs).tcell(ActiveIdx(std::get<0>(incident_nodes[n])));
      const mesh::CellGeometry<Cart2D::GDIM> cell_coords = tcell_view.coordinates();

      const auto node_coord = cell_coords.const_node_view(std::get<1>(incident_nodes[n]));

      const Real dx0  = node_coord[X0] - node_coord_ref[X0];
      const Real dx1  = node_coord[X1] - node_coord_ref[X1];
      const Real dist = std::sqrt(dx0 * dx0 + dx1 * dx1);

      BOOST_CHECK_CLOSE(dist, 0.0, 1.e-10);
    }
  }

  // Check that the arrays node_to_cell and cell_to_p1_node are consistent
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
          node_to_cell.const_block(p1_vert_pos[local_p1_node_idx]);
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
