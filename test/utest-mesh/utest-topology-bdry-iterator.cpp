/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bdry_topology_iterator_test
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

typedef Tria<Cart2D> MeshType;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(bdry_topology_iterator_validity_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

#if 0
  using bdry_topo_iter_t = MeshBoundarySet<Cart2D>::const_bdry_topo_iterator;

  const MeshBoundarySet<Cart2D> &mesh_bdry = mesh2d.all_boundaries();

  common::IteratorRange<bdry_topo_iter_t> bdry_range = mesh_bdry.domain_range("Bottom");

  bdry_topo_iter_t it1 = bdry_range.begin();
  bdry_topo_iter_t it2 = bdry_range.begin();

  bool result = (it1 == it2);
  // std::cout << "Result = " << result << " [should be true]" << std::endl;
  BOOST_CHECK_EQUAL(result, true);

  ++it2;
  result = (it1 == it2);
  BOOST_CHECK_EQUAL(result, false);
  // std::cout << "Result = " << result << " [should be false]" << std::endl;

  --it2;

  const bdry_topo_iter_t::value_type cell1 = *it1;
  const bdry_topo_iter_t::value_type cell2 = *(it2.operator->());

  std::cout << "Cell 1 = " << cell1 << std::endl;
  std::cout << "Cell 2 = " << cell2 << std::endl;
#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_topology_iterator_loop_utest)
{
  // std::cout << std::endl << "************ RUNNING CELL ITERATOR TEST
  // ************" << std::endl;

  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

#if 0
  /*
  std::cout << std::endl;
  std::cout
      << "*****************************************************************"
      << std::endl;
  std::cout << "Iterating over the boundary topology:" << std::endl;
  std::cout
      << "*****************************************************************"
      << std::endl;

  using bdry_topo_iter_t = MeshBoundarySet<Cart2D>::const_bdry_topo_iterator;
  // const MeshBoundarySet<Cart2D> &mesh_bdry = mesh2d.all_boundaries();

  common::IteratorRange<bdry_topo_iter_t> bdry_range =
      mesh2d.all_boundaries().domain_range("Bottom");

  for (bdry_topo_iter_t it = bdry_range.begin(); it != bdry_range.end(); ++it) {
    const bdry_topo_iter_t::value_type bdry_cell = *it;
    std::cout << bdry_cell << std::endl;
    std::cout << bdry_cell.coordinates() << std::endl;
  }
  */
#endif
}

// ----------------------------------------------------------------------------
