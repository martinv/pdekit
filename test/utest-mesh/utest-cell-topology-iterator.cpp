/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_topology_iterator_test
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

BOOST_AUTO_TEST_CASE(cell_topology_iterator_validity_utest)
{
  MeshType mesh2d("mesh2d");

  const std::string infilename = "rectangle_mixed_elem_p1.msh";

  gmsh::GmshReader meshreader;
  meshreader.read_mesh_from_file(infilename, mesh2d, "geo_dofs");

  typedef MeshType::const_cell_iterator cell_iterator_type;

  cell_iterator_type it1 = mesh2d.begin_cells();
  cell_iterator_type it2 = mesh2d.begin_cells();

  bool result = (it1 == it2);
  // std::cout << "Result = " << result << " [should be true]" << std::endl;
  BOOST_CHECK_EQUAL(result, true);

  ++it2;
  result = (it1 == it2);
  BOOST_CHECK_EQUAL(result, false);
  // std::cout << "Result = " << result << " [should be false]" << std::endl;

  --it2;

  const CellTopologyView<Cart2D> cell1 = *it1;
  const CellTopologyView<Cart2D> cell2 = *(it2.operator->());

  std::cout << "Cell 1 = " << cell1 << std::endl;
  std::cout << "Cell 2 = " << cell2 << std::endl;
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

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cell topology:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  for (MeshType::const_cell_iterator it = mesh2d.begin_cells(); it != mesh2d.end_cells(); ++it)
  {
    const CellTopologyView<Cart2D> topo_cell = *it;
    std::cout << topo_cell << std::endl;
  }

  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cell topology, only active cells:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  for (MeshType::const_active_cell_iterator it = mesh2d.begin_cells_active();
       it != mesh2d.end_cells_active(); ++it)
  {
    const CellTopologyView<Cart2D> topo_cell = *it;
    std::cout << topo_cell << std::endl;
  }
}

// ================================ TIMINGS ===================================

BOOST_AUTO_TEST_CASE(cell_topology_iterator_timings_utest)
{
  MeshType mesh2d("mesh2d");

  mesh::MeshCreator::make_unit_quad(mesh2d, "geo_dofs", 1001,
                                    false); /// 1 million quads

  std::cout << "Number of cells in mesh = " << mesh2d.nb_active_cells() << std::endl;

  clock_t start, end;
  Real elapsed;

  // Loop over cells based on their types
  std::cout << std::endl;
  std::cout << "*****************************************************************" << std::endl;
  std::cout << "Iterating over the cell topology:" << std::endl;
  std::cout << "*****************************************************************" << std::endl;

  start = clock();
  for (MeshType::const_cell_iterator it = mesh2d.begin_cells(); it != mesh2d.end_cells(); ++it)
  {
    const CellTopologyView<Cart2D> topo_cell             = *it;
    const common::ArrayView<const Uint, _1D, Uint> faces = topo_cell.incident_facets();
    const StdRegion cell_type                            = topo_cell.std_region();
    // BOOST_CHECK_EQUAL(faces.size(), 4);
    // std::cout << topo_cell << std::endl;
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout << "1 loop over all topology cells: " << elapsed << " sec" << std::endl;
}

// ----------------------------------------------------------------------------
