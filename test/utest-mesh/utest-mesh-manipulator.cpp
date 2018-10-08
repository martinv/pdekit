/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_manipulator_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "mesh/Tria.hpp"
#include "mesh/io/MeshManipulator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType2D;
typedef Tria<Cart3D> MeshType3D;

struct MeshManipulatorFixture
{
  gmsh::GmshReader meshreader;
  gmsh::GmshWriter meshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MeshManipulator_TestSuite, MeshManipulatorFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(scale_mesh_2D)
{
  MeshType2D mesh2d("Mesh2D");

  // ----------------------------------------------------------------------------

  const std::string infilename2d = "test_p1_tri.msh";
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  const Real scale_factor = 2.0;
  MeshManipulator::scale_mesh(mesh2d, "geo_dofs", scale_factor);

  const std::string outfilename2d = "scaled_" + infilename2d;
  meshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(translate_mesh_2D)
{
  MeshType2D mesh2d("Mesh2D");

  // ----------------------------------------------------------------------------

  const std::string infilename2d = "test_p1_tri.msh";
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  math::DenseSVec<Real, Cart2D::GDIM> trans_vec;
  trans_vec[X0] = 0.1;
  trans_vec[X1] = -0.2;

  MeshManipulator::translate_mesh(mesh2d, "geo_dofs", trans_vec);

  const std::string outfilename2d = "translated_" + infilename2d;
  meshwriter.write_mesh_to_file(mesh2d, "geo_dofs", outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
