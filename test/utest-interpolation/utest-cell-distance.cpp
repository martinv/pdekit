/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE cell_distance_utest
#include <boost/test/unit_test.hpp>
#include <memory>

#include "interpolation/mesh_function/MeshFunctionTools.hpp"

#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;

typedef mesh::Cart2D MeshConfig2D;
typedef mesh::Cart3D MeshConfig3D;

typedef mesh::Tria<MeshConfig2D> MeshType2D;
typedef mesh::Tria<MeshConfig3D> MeshType3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_distance_2D_utest)
{
  MeshType2D::shared_ptr mesh2d = std::make_shared<MeshType2D>("square2D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_square_tri_p2.msh", *mesh2d, "geo_dofs");

  interpolation::ScalarMeshFunction<Uint> distance("", "distance");
  std::vector<mesh::ActiveIdx> seed_cells = {mesh::ActiveIdx(200), mesh::ActiveIdx(300),
                                             mesh::ActiveIdx(100)};
  common::ArrayView<const mesh::ActiveIdx, _1D, Uint> seeds(seed_cells.data(), seed_cells.size());

  interpolation::MeshFunctionTools::compute_cell_distance(*mesh2d, seeds, distance);

  mesh::gmsh::GmshWriter mesh_writer;
  mesh_writer.write_mesh_to_file(*mesh2d, "geo_dofs", "distance_2D.msh");
  mesh_writer.append_cell_function_to_file(*mesh2d, "distance_2D.msh", distance, "distance");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_distance_3D_utest)
{
  MeshType3D::shared_ptr mesh3d = std::make_shared<MeshType3D>("mesh3D");

  mesh::gmsh::GmshReader mesh_reader;
  mesh_reader.read_mesh_from_file("unit_cube_tet_p3.msh", *mesh3d, "geo_dofs");

  interpolation::ScalarMeshFunction<Uint> distance("", "distance");
  std::vector<mesh::ActiveIdx> seed_cells = {mesh::ActiveIdx(7000)};
  common::ArrayView<const mesh::ActiveIdx, _1D, Uint> seeds(seed_cells.data(), seed_cells.size());

  interpolation::MeshFunctionTools::compute_cell_distance(*mesh3d, seeds, distance);

  mesh::gmsh::GmshWriter mesh_writer;
  mesh_writer.write_mesh_to_file(*mesh3d, "geo_dofs", "distance_3D.msh");
  mesh_writer.append_cell_function_to_file(*mesh3d, "distance_3D.msh", distance, "distance");
}

// ----------------------------------------------------------------------------
