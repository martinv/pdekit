/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_path_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <iostream>

/// PDEKIT headers
#include "mesh/MeshConstants.hpp"
#include "mesh/containers/CellPath.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(path_segment_utest)
{
  detail::PathSegment segment;

  segment.set_entry(0, 4);
  segment.set_entry(1, 2);
  segment.set_entry(2, 1);
  segment.set_entry(3, 3);
  segment.set_entry(4, 4);
  segment.set_entry(5, 3);
  segment.set_entry(6, 0);
  segment.set_entry(7, 4);

  BOOST_CHECK_EQUAL(segment.entry(0), 4u);
  BOOST_CHECK_EQUAL(segment.entry(1), 2u);
  BOOST_CHECK_EQUAL(segment.entry(2), 1u);
  BOOST_CHECK_EQUAL(segment.entry(3), 3u);
  BOOST_CHECK_EQUAL(segment.entry(4), 4u);
  BOOST_CHECK_EQUAL(segment.entry(5), 3u);
  BOOST_CHECK_EQUAL(segment.entry(6), 0u);
  BOOST_CHECK_EQUAL(segment.entry(7), 4u);

  segment.set_entry(3, 0);
  BOOST_CHECK_EQUAL(segment.entry(3), 0u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_path_construction_utest)
{
  mesh::CellPath cp1;
  BOOST_CHECK_EQUAL(cp1.base_cell_id(), mesh::FlatIdx(INVALID_CELL_ID));
  BOOST_CHECK_EQUAL(cp1.capacity(), 0u);
  BOOST_CHECK_EQUAL(cp1.size(), 0u);

  const std::vector<Uint> path_entries = {2, 1, 1, 3, 0, 2, 0, 4, 2, 9};
  mesh::CellPath cp2                   = mesh::CellPath(mesh::FlatIdx(150), path_entries);

  BOOST_CHECK_EQUAL(cp2.base_cell_id(), mesh::FlatIdx(150));
  BOOST_CHECK_EQUAL(cp2.capacity(), 15u);
  BOOST_CHECK_EQUAL(cp2.size(), 10u);

  for (Uint i = 0; i < path_entries.size(); ++i)
  {
    BOOST_CHECK_EQUAL(cp2.path_entry(i), path_entries[i]);
  }

  std::cout << cp2 << std::endl;

  mesh::CellPath cp3(mesh::FlatIdx(100), {});
  BOOST_CHECK_EQUAL(cp3.base_cell_id(), mesh::FlatIdx(100));
  BOOST_CHECK_EQUAL(cp3.capacity(), 0u);
  BOOST_CHECK_EQUAL(cp3.size(), 0u);

  std::cout << cp3 << std::endl;

  CellPath cp4;
  cp4 = cp2;
  std::cout << "cp2" << std::endl;
  std::cout << cp2 << std::endl;
  std::cout << "Capacity of cp2 = " << cp2.capacity() << std::endl;
  std::cout << "cp4" << std::endl;
  std::cout << cp4 << std::endl;
  cp2.print();
  cp4.print();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_path_filling_utest)
{
  mesh::CellPath path;

  // Pass 1

  for (Uint i = 0; i < detail::PathSegment::max_num_entries - 1; ++i)
  {
    path.push_back(i + 1);
  }

  std::cout << path << std::endl;

  for (Uint i = 0; i < detail::PathSegment::max_num_entries - 1; ++i)
  {
    BOOST_CHECK_EQUAL(path.path_entry(i), i + 1);
  }

  BOOST_CHECK_EQUAL(path.size(), detail::PathSegment::max_num_entries - 1);
  BOOST_CHECK_EQUAL(path.capacity(), detail::PathSegment::max_num_entries - 1);

  // Pass 2

  for (Uint i = 0; i < detail::PathSegment::max_num_entries; ++i)
  {
    path.push_back(detail::PathSegment::max_entry_value - 1 - i);
  }

  std::cout << path << std::endl;

  for (Uint i = 0; i < detail::PathSegment::max_num_entries; ++i)
  {
    BOOST_CHECK_EQUAL(path.path_entry(detail::PathSegment::max_num_entries - 1 + i),
                      detail::PathSegment::max_entry_value - 1 - i);
  }

  BOOST_CHECK_EQUAL(path.size(), 2 * detail::PathSegment::max_num_entries - 1);
  BOOST_CHECK_EQUAL(path.capacity(), 2 * detail::PathSegment::max_num_entries - 1);

  // Pass 3

  for (Uint i = 0; i < detail::PathSegment::max_num_entries; ++i)
  {
    path.push_back(2 + i);
  }

  std::cout << path << std::endl;

  for (Uint i = 0; i < detail::PathSegment::max_num_entries - 1; ++i)
  {
    BOOST_CHECK_EQUAL(path.path_entry(2 * detail::PathSegment::max_num_entries - 1 + i), 2 + i);
  }

  BOOST_CHECK_EQUAL(path.size(), 3 * detail::PathSegment::max_num_entries - 1);
  BOOST_CHECK_EQUAL(path.capacity(), 3 * detail::PathSegment::max_num_entries - 1);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(cell_path_emptying_utest)
{
  mesh::CellPath path;

  // Pass 1

  for (Uint i = 0; i < detail::PathSegment::max_num_entries - 1; ++i)
  {
    path.push_back(i + 1);
  }

  path.push_back(11);

  std::cout << path << std::endl;
  std::cout << "\nStarting removal" << std::endl;
  for (Uint i = 0; i < detail::PathSegment::max_num_entries - 1; ++i)
  {
    path.pop();
    std::cout << path << std::endl;
  }
  path.pop();
  std::cout << path << std::endl;

  path.push_back(15);
  path.push_back(14);
  std::cout << path << std::endl;
}

// ----------------------------------------------------------------------------
