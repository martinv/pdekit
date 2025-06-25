/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE vtk_reader_writer_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"
#include "mesh/io/vtk/VtkReader.hpp"
#include "mesh/io/vtk/VtkWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> MeshType2D;
typedef Tria<Cart3D> MeshType3D;
typedef Tria<Surf2D3D> MeshTypeSurf3D;

struct VtkReaderWriterFixture
{
  vtk::VtkReader vtk_reader;
  vtk::VtkWriter vtk_writer;
  gmsh::GmshWriter gmsh_writer;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(VtkReaderWriterAdaptive_TestSuite, VtkReaderWriterFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vtk_reader_writer_cart_3D)
{
  MeshType3D mesh3d("Mesh3D");

  std::cout << std::endl << "************ RUNNING VTK Reader TEST (3D) ************" << std::endl;

  const std::string infilename3d = "Cube.vtk";
  vtk_reader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  common::PtrHandle<mesh::DofMap<Cart3D>> geo_dofs = mesh3d.dof_storage("geo_dofs");

  using skeleton_iter_t        = mesh::Tria<Cart3D>::const_trace_dof_iterator;
  skeleton_iter_t facets_begin = mesh3d.cbegin_skeleton_dofs(_2D, *geo_dofs);
  skeleton_iter_t facets_end   = mesh3d.cend_skeleton_dofs(_2D, *geo_dofs);

  std::ofstream outfile;
  outfile.setf(std::ios::scientific);

  outfile.open("skin.stl");
  outfile << "solid surface\n";

  for (skeleton_iter_t facet_it = facets_begin; facet_it != facets_end; ++facet_it)
  {
    if (facet_it->size() == 1)
    {
      const mesh::CellGeometry<_3D> facet_geo          = facet_it->cell_geometry(0);
      const mesh::CellGeometry<_3D>::node_view_t node0 = facet_geo.const_node_view(0);
      const mesh::CellGeometry<_3D>::node_view_t node1 = facet_geo.const_node_view(1);
      const mesh::CellGeometry<_3D>::node_view_t node2 = facet_geo.const_node_view(2);

      const Real v0x = node0[X0];
      const Real v0y = node0[X1];
      const Real v0z = node0[X2];

      const Real v1x = node1[X0];
      const Real v1y = node1[X1];
      const Real v1z = node1[X2];

      const Real v2x = node2[X0];
      const Real v2y = node2[X1];
      const Real v2z = node2[X2];

      const Real e0x = v1x - v0x;
      const Real e0y = v1y - v0y;
      const Real e0z = v1z - v0z;

      const Real e1x = v2x - v0x;
      const Real e1y = v2y - v0y;
      const Real e1z = v2z - v0z;

      const Real nx = e0y * e1z - e1y * e0z;
      const Real ny = e1x * e0z - e0x * e1z;
      const Real nz = e0x * e1y - e1x * e0y;

      outfile << "facet normal " << nx << " " << ny << " " << nz << "\n";
      outfile << "    outer loop\n";
      outfile << "        vertex " << v0x << " " << v0y << " " << v0z << "\n";
      outfile << "        vertex " << v1x << " " << v1y << " " << v1z << "\n";
      outfile << "        vertex " << v2x << " " << v2y << " " << v2z << "\n";
      outfile << "    endloop\n";
      outfile << "endfacet\n";
    }
  }

  outfile << "endsolid";
  outfile.close();

  const std::string outfilename3d = "output_intercostal.msh";
  gmsh_writer.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d);

  std::cout << std::endl << "************ FINISHED VTK Reader TEST (3D) ************" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vtk_reader_writer_surface_2D3D)
{
  MeshTypeSurf3D mesh3d("Mesh3D");

  std::cout << std::endl
            << "************ RUNNING VTK Reader TEST (manifold 3D) ************" << std::endl;

  const std::string infilename3d_polydata = "Leaflet_polygon.vtk";
  vtk_reader.read_mesh_from_file(infilename3d_polydata, mesh3d, "geo_dofs");

  const std::string outfilename3d_polydata = "output_leaflet_polygon.msh";
  gmsh_writer.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d_polydata);

  std::cout << std::endl
            << "************ FINISHED VTK Reader TEST (manifold 3D) ************" << std::endl;
}

// ----------------------------------------------------------------------------

/*
BOOST_AUTO_TEST_CASE(gmsh_reader_writer_3D)
{

  MeshType3D mesh3d("Mesh3D");

  //
----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Reader TEST (3D)
************" << std::endl;

  const std::string infilename3d = "bump_p1_tet.msh";
  meshreader.read_mesh_from_file(infilename3d, mesh3d, "geo_dofs");

  std::cout << std::endl << "************ FINISHED GMSH Reader TEST (3D)
************" << std::endl;

  //
----------------------------------------------------------------------------

  std::cout << std::endl << "************ RUNNING GMSH Writer TEST (3D)
************" << std::endl;

  const std::string outfilename3d = "output_" + infilename3d;
  meshwriter.write_mesh_to_file(mesh3d, "geo_dofs", outfilename3d);

  std::cout << std::endl << "************ FINISHED GMSH Writer TEST (3D)
************" << std::endl;
}
*/
// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
