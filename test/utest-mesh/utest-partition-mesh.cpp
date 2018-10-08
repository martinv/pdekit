/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_partitioning_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "graph/GraphPartitioner.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType2D;
typedef Tria<Cart3D> MeshType3D;

struct MeshParitioningFixture
{
  gmsh::GmshReader meshreader;
  gmsh::GmshWriter meshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MeshPartitioning_TestSuite, MeshParitioningFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(partion_mesh_2D_scotch)
{
  MeshType2D mesh2d("Mesh2D");

  // ----------------------------------------------------------------------------

  const std::string infilename2d = "test_p3_tri.msh";
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  common::BlockArray<Uint, Uint> mesh_graph;
  std::vector<Uint> partition_ids;

  mesh2d.build_dual_graph_undirected(mesh_graph);

  graph::GraphPartitioner partitioner;

  const Uint nb_parts = 5;

  partitioner.part_graph(mesh_graph, nb_parts, "scotch", partition_ids);

  const common::ArrayView<const Uint, _1D, Uint> mesh_part_ids(partition_ids.data(),
                                                               partition_ids.size());

  const std::string outfilename2d = "partitioned_" + infilename2d;

  meshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
  meshwriter.append_cell_function_to_file(mesh2d, outfilename2d, mesh_part_ids, "part_id");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(gmsh_reader_writer_3D)
{

  MeshType3D mesh3d("Mesh3D");

  const std::string infilename3d = "bump_p1_tet.msh";
  meshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  common::BlockArray<Uint, Uint> mesh_graph;
  std::vector<Uint> partition_ids;

  mesh3d.build_dual_graph_undirected(mesh_graph);

  graph::GraphPartitioner partitioner;

  const Uint nb_parts = 5;

  partitioner.part_graph(mesh_graph, nb_parts, "scotch", partition_ids);

  const common::ArrayView<const Uint, _1D, Uint> mesh_part_ids(partition_ids.data(),
                                                               partition_ids.size());

  const std::string outfilename3d = "partitioned_" + infilename3d;

  meshwriter.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d);
  meshwriter.append_cell_function_to_file(mesh3d, outfilename3d, mesh_part_ids, "part_id");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
