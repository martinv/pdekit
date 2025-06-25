/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE skeleton_iterator_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "mesh/io/MeshCreator.hpp"
#include "mesh/io/gmsh/GmshReader.hpp"
#include "mesh/io/gmsh/GmshWriter.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

using MeshType = Tria<Cart2D>;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(skeleton_iterator_validity_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

  using skeleton_iterator_t = MeshType::const_skeleton_iterator;

  skeleton_iterator_t it1 = mesh2d.cbegin_skeleton(_1D);
  skeleton_iterator_t it2 = mesh2d.cbegin_skeleton(_1D);

  bool result = (it1 == it2);
  // std::cout << "Result = " << result << " [should be true]" << std::endl;
  BOOST_CHECK_EQUAL(result, true);

  ++it2;
  result = (it1 == it2);
  BOOST_CHECK_EQUAL(result, false);

  // std::cout << "Result = " << result << " [should be false]" << std::endl;

  --it2;

  const TraceIncidences li1 = it1->incidences();
  const TraceIncidences li2 = (*(it2.operator->())).incidences();

  std::cout << "Local incidences 1 = " << li1 << std::endl;
  std::cout << "Local incidences 2 = " << li2 << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(skeleton_iterator_loop_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename2d = "unit_square_mini.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");

  using skeleton_iterator_t = MeshType::const_skeleton_iterator;

  for (skeleton_iterator_t it = mesh2d.cbegin_skeleton(_1D); it != mesh2d.cend_skeleton(_1D); ++it)
  {
    std::cout << it->incidences() << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dof_skeleton_iterator_loop_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename2d = "unit_square_mini.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename2d, mesh2d, "geo_dofs");
  common::PtrHandle<MeshType::dof_storage_type> geo_dofs = mesh2d.dof_storage("geo_dofs");
  (*geo_dofs).upgrade(mesh2d, P2, PointSetID::Equidist);

  using trace_dof_iterator_t = MeshType::const_trace_dof_iterator;

  for (trace_dof_iterator_t it = mesh2d.cbegin_skeleton_dofs(_1D, *geo_dofs);
       it != mesh2d.cend_skeleton_dofs(_1D, *geo_dofs); ++it)
  {
    std::cout << "--------------------------------" << std::endl;
    const auto view = *it;
    for (Uint i = 0; i < view.size(); ++i)
    {
      std::cout << "Topology cell type: " << view.geo_cell_type(i).as_string() << std::endl;
      std::cout << "Dof cell type: " << view.cell_type(i).as_string() << std::endl;
      std::cout << view.mesh_entity(i) << std::endl;
      std::cout << "Topology face geometry:" << std::endl;
      std::cout << view.cell_geometry(i) << std::endl;
    }
  }
}

// ----------------------------------------------------------------------------
