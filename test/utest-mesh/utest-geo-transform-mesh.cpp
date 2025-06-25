/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_geometry_transformation_test
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

struct GeoTransformMeshFixture
{
  gmsh::GmshReader gmshreader;
  gmsh::GmshWriter gmshwriter;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(GeoTransformMesh_TestSuite, GeoTransformMeshFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(geo_transform_mesh_2D)
{
  std::cout << "\n************ RUNNING MESH GEOMETRY TRANSFORM TEST - 2D "
               "************"
            << std::endl;

  // Create a mesh

  Tria<Cart2D> mesh2d("mesh2d");

  // Read some data in it
  const std::string infilename2d = "test_p1_tri.msh";
  gmshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  mesh2d.geo_transform(
      [](const math::DenseConstVecView<Real> &coord_in, math::DenseVecView<Real> &coord_out) {
        coord_out[X0] = coord_in[X0] + 0.1;
        coord_out[X1] = coord_in[X1] + 1.0;
      });

  // Write the mesh in the output file
  const std::string outfilename2d = "mesh_geo_transformed_" + infilename2d;
  gmshwriter.write_mesh_to_file(mesh2d, outfilename2d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(geo_transform_mesh_3D)
{
  std::cout << "\n************ RUNNING MESH GEOMETRY TRANSFORM TEST - 3D "
               "************"
            << std::endl;

  Tria<Cart3D> mesh3d("mesh3d");

  // Read some data in it
  const std::string infilename3d = "bump_p1_tet.msh";
  gmshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  mesh3d.geo_transform(
      [](const math::DenseConstVecView<Real> &coord_in, math::DenseVecView<Real> &coord_out) {
        coord_out[X0] = coord_in[X0] + 0.1;
        coord_out[X1] = coord_in[X1] + 1.0;
        coord_out[X2] = coord_in[X2] - 1.3;
      });

  // Write the mesh in the output file
  const std::string outfilename3d = "mesh_geo_transformed_" + infilename3d;
  gmshwriter.write_mesh_to_file(mesh3d, outfilename3d);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
