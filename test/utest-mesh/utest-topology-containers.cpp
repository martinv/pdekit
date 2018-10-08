/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE adaptive_mesh_unit_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers

#include "mesh/MeshConfig.hpp"
#include "mesh/containers/DofMap.hpp"
#include "mesh/containers/TriaCells.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_connectivity_utest)
{
  internal::TriaCells<Cart2D> cell_connectivity;

  // Faces of 6 triangles and one qudrilateral: total = 6 x 3 + 1 x 4 = 22
  // faces
  std::vector<Uint> cell_faces = {1, 3, 2,  2,  8,  7,  4,  5, 6,  6,  9,
                                  3, 8, 12, 11, 14, 10, 11, 9, 13, 15, 12};
  std::vector<StdRegion> cell_types;

  BOOST_CHECK_EQUAL(cell_faces.size(), 22U);

  StdRegion stdr;

  for (Uint i = 0; i < 6; ++i)
  {
    stdr.change_type(ElemShape::Triag, P1, PointSetID::Equidist);
    cell_types.push_back(stdr);
  }

  stdr.change_type(ElemShape::Quad, P1, PointSetID::Equidist);
  cell_types.push_back(stdr);

  common::BlockArray<Real, Uint> cell_coordinates;
  cell_coordinates.resize(6 * 3 * 2 + 4 * 2, 7);

  cell_connectivity.emplace_cells(cell_types, cell_faces, cell_coordinates);

  std::cout << cell_connectivity << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(dof_handler_utest)
{
  std::cout << "Creating dof handler" << std::endl;

  DofMap<Cart2D> dof_handler(nullptr, "my_dof_handler");
  dof_handler.init(_2D);

  std::vector<Uint> dof_ids = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::vector<IncidenceEntry> cell_ids;
  std::vector<PointSetTag> type_ids;
  for (Uint i = 0; i < 3; ++i)
  {
    cell_ids.push_back(IncidenceEntry(i, 0u));
    type_ids.push_back(PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist));
  }

  /*
  dof_handler.emplace_dof_data(dof_ids, cell_ids, type_ids);
  dof_handler.print_cell_types();

  send_to_output_stream(std::cout, dof_handler);
  dof_handler.print_cell_types();
  std::cout << std::endl;
  */
}

// ----------------------------------------------------------------------------
