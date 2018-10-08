/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_upgrade_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/Constants.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

struct UpgradeMeshUtestFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(UpgradeMesh_TestSuite, UpgradeMeshUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(upgrade_mesh_2D)
{
  std::cout << "\n************ RUNNING MESH UPGRADE TEST - 2D ************" << std::endl;

  // Create a mesh

  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  mesh2d.change_std_region_types([](const mesh::PointSetTag &tag_in) {
    return mesh::PointSetTag(tag_in.elem_shape(), tag_in.poly_order() + 2, tag_in.ref_topology());
  });

  // Write the mesh in the output file
  const std::string outfilename2d = "mesh_upgraded_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(upgrade_mesh_3D)
{
  std::cout << "\n************ RUNNING MESH UPGRADE TEST - 3D ************" << std::endl;

  Tria<Cart3D> mesh3d("mesh3d");

  // Read some data in it
  const std::string infilename3d = "bump_p1_tet.msh";
  gmshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  mesh3d.change_std_region_types([](const mesh::PointSetTag &tag_in) {
    return mesh::PointSetTag(tag_in.elem_shape(), tag_in.poly_order() + 2, tag_in.ref_topology());
  });

  // Write the mesh in the output file
  const std::string outfilename3d = "mesh_upgraded_" + infilename3d;
  gmshwriter.write_mesh_to_file(mesh3d, outfilename3d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(upgrade_dofs_2D)
{
  std::cout << "\n************ RUNNING DOF UPGRADE TEST - 2D ************" << std::endl;

  // Create a mesh

  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri.msh";

  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");
  common::PtrHandle<Tria<Cart2D>::dof_storage_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  (*geo_dofs).upgrade(mesh2d, P5, PointSetID::Warpblend);

  // Write the mesh in the output file
  const std::string outfilename2d = "dof_upgraded_" + infilename2d;

  gmshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(upgrade_dofs_3D)
{
  std::cout << "\n************ RUNNING DOF UPGRADE TEST - 3D ************" << std::endl;

  Tria<Cart3D> mesh3d("mesh3d");

  // Read some data in it
  const std::string infilename3d = "bump_p1_tet.msh";

  gmshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  // Perform the actual upgrade
  common::PtrHandle<Tria<Cart3D>::dof_storage_type> geo_dofs = mesh3d.dof_storage("geo_dofs");
  (*geo_dofs).upgrade(mesh3d, P4, PointSetID::Equidist);

  // Write the mesh in the output file
  const std::string outfilename3d = "dof_upgraded_" + infilename3d;

  gmshwriter.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(clone_dofs_2D)
{
  std::cout << "\n************ RUNNING TOPOLOGY CLONING TEST - 2D ************" << std::endl;

  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri.msh";

  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "continuous_dofs");
  common::PtrHandle<Tria<Cart2D>::dof_storage_type> continuous_dofs =
      mesh2d.dof_storage("continuous_dofs");
  common::PtrHandle<Tria<Cart2D>::dof_storage_type> discontinuous_dofs =
      mesh2d.create_dof_storage("discontinuous_dofs");

  // Perform the actual upgrade
  Tria<Cart2D>::dof_storage_type::clone_continuous(mesh2d, *continuous_dofs, *discontinuous_dofs,
                                                   P2, PointSetID::Equidist);

  // Tria<Cart2D>::clone(mesh2d_in, mesh2d_out, 4, Warpblend);
  const Uint nb_edges = mesh2d.active_skeleton_size(_1D);

  BOOST_CHECK_EQUAL(nb_edges, 80u);

  // Write the mesh in the output file
  const std::string outfilename2d = "cloned_" + infilename2d;

  gmshwriter.write_mesh_to_file(mesh2d, "discontinuous_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(clone_dofs_3D)
{
  std::cout << "\n************ RUNNING TOPOLOGY CLONING TEST - 3D ************" << std::endl;

  Tria<Cart3D> mesh3d("mesh3d");

  // Read some data in it
  const std::string infilename3d = "bump_p1_tet.msh";

  gmshreader.read_mesh_from_file(infilename3d, mesh3d, "continuous_dofs");
  common::PtrHandle<Tria<Cart3D>::dof_storage_type> continuous_dofs =
      mesh3d.dof_storage("continuous_dofs");
  common::PtrHandle<Tria<Cart3D>::dof_storage_type> discontinuous_dofs =
      mesh3d.create_dof_storage("discontinuous_dofs");

  // Perform the actual upgrade
  Tria<Cart3D>::dof_storage_type::clone_continuous(mesh3d, *continuous_dofs, *discontinuous_dofs,
                                                   P3, PointSetID::Equidist);

  const Uint nb_nodes = (*continuous_dofs).nb_nodes();
  const Uint nb_edges = mesh3d.active_skeleton_size(_1D);
  const Uint nb_faces = mesh3d.active_skeleton_size(_2D);
  const Uint nb_cells = mesh3d.nb_active_cells();

  BOOST_CHECK_EQUAL(nb_nodes, 2317u);
  BOOST_CHECK_EQUAL(nb_edges, 13068u);
  BOOST_CHECK_EQUAL(nb_faces, 20080u);
  BOOST_CHECK_EQUAL(nb_cells, 9328u);

  const Uint LHS = nb_cells + nb_edges + 1;
  const Uint RHS = nb_nodes + nb_faces;
  BOOST_CHECK_EQUAL(LHS, RHS);

  // Write the mesh in the output file
  const std::string outfilename3d = "cloned_" + infilename3d;

  gmshwriter.write_mesh_to_file(mesh3d, "discontinuous_dofs", outfilename3d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
