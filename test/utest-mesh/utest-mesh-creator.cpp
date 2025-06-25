/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_creator_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
// using namespace pdekit::common;

struct MeshCreatorFixture
{
  /// common setup for each test case
  MeshCreatorFixture()
  {
  }

  ~MeshCreatorFixture()
  {
  }

  // MeshCreator m_creator;
  gmsh::GmshWriter m_gmsh_writer;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MeshCreator_TestSuite, MeshCreatorFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_creator_triangle_2D)
{
  Tria<Cart2D> unit_triag_mesh("unit_triag_mesh");
  MeshCreator::make_unit_triangle(unit_triag_mesh, "geo_dofs", 50);

  m_gmsh_writer.write_mesh_to_file(unit_triag_mesh, "geo_dofs", "created_triag_mesh.msh");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_creator_quad_2D)
{
  Tria<Cart2D> unit_quad_mesh("unit_quad_mesh");
  MeshCreator::make_unit_quad(unit_quad_mesh, "geo_dofs", 50, false);

  m_gmsh_writer.write_mesh_to_file(unit_quad_mesh, "geo_dofs", "created_quad_mesh.msh");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
