/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_geometry_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>
#include <memory>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "mesh/DofCoordinates.hpp"
#include "mesh/MeshPredicates.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

typedef Tria<Cart2D> Mesh2D;
typedef Tria<Cart3D> Mesh3D;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mesh_geometry_utest)
{
  std::cout << std::endl << "************ RUNNING MESH GEOMETRY TEST 2D ************" << std::endl;

  Mesh2D::shared_ptr mesh2d = std::make_shared<Mesh2D>("mesh2d");

  // Mesh root("root");

  std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, *mesh2d, "geo_dofs_2D");

  // const Mesh2D::cell_connectivity& connectivity = root->topology().cells();

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cells group by group:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  const Mesh2D::dof_storage_type &cells2d = *mesh2d->dof_storage("geo_dofs_2D");

  for (Mesh2D::dof_storage_type::const_dof_range_typed cells : cells2d.all_active_dof_groups())
  {
    for (Mesh2D::dof_storage_type::const_dof_iterator_typed cell_iter = cells.begin();
         cell_iter != cells.end(); ++cell_iter)
    {
      // std::cout << *cell_iter << std::endl;

      const CellTopologyView<Cart2D> tview = cell_iter->tcell();
      const MeshEntity cell                = cell_iter->mesh_entity();
      CellGeometry<Cart2D::GDIM> cc        = tview.coordinates();

      std::cout << "[" << cell.idx() << "]  " << cell << std::endl;

      for (Uint n = 0; n < cc.size(); ++n)
      {
        std::cout << "\t" << cc.const_node_view(n) << std::endl;
      }

      std::cout << "\t\tSub-element coordinates:" << std::endl;
      for (Uint i_sub = 0; i_sub < cell.nb_sub_elements(_1D); ++i_sub)
      {
        const mesh::CellGeometry<Cart2D::GDIM> sub_coord = tview.coordinates(_1D, i_sub);
        for (Uint n = 0; n < sub_coord.size(); ++n)
        {
          std::cout << "\t\t" << sub_coord.const_node_view(n) << std::endl;
        }
        std::cout << std::endl;
      }

      std::cout << std::endl;
    }
  }

  std::cout << std::endl << "************ FINISHED MESH GEOMETRY TEST 2D ************" << std::endl;

  std::cout << std::endl << "************ RUNNING MESH GEOMETRY TEST 3D ************" << std::endl;

  Mesh3D::shared_ptr mesh3d = std::make_shared<Mesh3D>("mesh3d");

  infilename = "cube.msh";

  meshreader.read_mesh_from_file(infilename, *mesh3d, "geo_dofs_3D");

  // const Mesh2D::cell_connectivity& connectivity = root->topology().cells();

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cells group by group:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  const Mesh3D::dof_storage_type &cells3d = *mesh3d->dof_storage("geo_dofs_3D");

  for (Mesh3D::dof_storage_type::const_dof_range_typed cells : cells3d.all_active_dof_groups())
  {
    for (Mesh3D::dof_storage_type::const_dof_iterator_typed cell_iter = cells.begin();
         cell_iter != cells.end(); ++cell_iter)
    {
      // std::cout << *cell_iter << std::endl;
      const CellTopologyView<Cart3D> tview = cell_iter->tcell();
      const MeshEntity cell                = cell_iter->mesh_entity();
      CellGeometry<Cart3D::GDIM> cc        = tview.coordinates();

      std::cout << "[" << cell.idx() << "]  " << cell << std::endl;

      for (Uint n = 0; n < cc.size(); ++n)
      {
        std::cout << "\t" << cc.const_node_view(n) << std::endl;
      }

      std::cout << "\t\tSub-element coordinates:" << std::endl;
      for (Uint i_sub = 0; i_sub < cell.nb_sub_elements(_2D); ++i_sub)
      {
        const mesh::CellGeometry<Cart3D::GDIM> sub_coord = tview.coordinates(_2D, i_sub);
        for (Uint n = 0; n < sub_coord.size(); ++n)
        {
          std::cout << "\t\t" << sub_coord.const_node_view(n) << std::endl;
        }
        std::cout << std::endl;
      }

      std::cout << std::endl;
    }
  }

  std::cout << std::endl << "************ FINISHED MESH GEOMETRY TEST 3D ************" << std::endl;
}

// ----------------------------------------------------------------------------
