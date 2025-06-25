/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gmsh_reader_writer_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/Tria.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType2D;
typedef Tria<Cart3D> MeshType3D;

struct GmshReaderWriterFixture
{
  gmsh::GmshReader meshreader;
  gmsh::GmshWriter meshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(GmshReaderWriterAdaptive_TestSuite, GmshReaderWriterFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(gmsh_reader_writer_2D)
{
  MeshType2D mesh2d("Mesh2D");

  // ----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Reader TEST (2D) ************" << std::endl;

  const std::string infilename2d = "test_p3_tri.msh";
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  std::cout << std::endl << "************ FINISHED GMSH Reader TEST (2D) ************" << std::endl;

  // ----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Writer TEST (2D) ************" << std::endl;

  const std::string outfilename2d = "output_" + infilename2d;
  meshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);

  std::cout << std::endl << "************ FINISHED GMSH Writer TEST (2D) ************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(gmsh_reader_writer_3D)
{

  MeshType3D mesh3d("Mesh3D");

  // ----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Reader TEST (3D) ************" << std::endl;

  const std::string infilename3d = "bump_p1_tet.msh";
  meshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  std::cout << std::endl << "************ FINISHED GMSH Reader TEST (3D) ************" << std::endl;

  // ----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Writer TEST (3D) ************" << std::endl;

  const std::string outfilename3d = "output_" + infilename3d;
  meshwriter.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d);

  std::cout << std::endl << "************ FINISHED GMSH Writer TEST (3D) ************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
