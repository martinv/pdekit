/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE graph_partitioner_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "graph/GraphPartitioner.hpp"

using namespace pdekit;

struct GraphPartitionerFixture
{
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(GraphPartitioner_TestSuite, GraphPartitionerFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(graph_partitioner)
{
  common::BlockArray<Uint, Uint> mesh_graph;

  std::unique_ptr<std::vector<Uint>> vert_ids(new std::vector<Uint>());
  std::unique_ptr<std::vector<Uint>> vert_adj_blocks(new std::vector<Uint>());

  (*vert_ids) = {1,  5,  0,  6,  2,  1,  7,  3,  2,  8,  4,  3,  9,  15, 0,  6,  10, 5,  11, 7,
                 1,  6,  12, 8,  2,  7,  13, 9,  3,  8,  14, 16, 4,  5,  11, 31, 10, 32, 12, 6,
                 11, 33, 13, 7,  12, 34, 14, 8,  13, 35, 17, 9,  4,  16, 15, 9,  17, 16, 14, 18,
                 17, 35, 19, 18, 25, 20, 19, 30, 22, 26, 31, 21, 27, 23, 32, 22, 28, 24, 33, 23,
                 29, 25, 34, 24, 30, 19, 35, 21, 27, 26, 22, 28, 27, 23, 29, 28, 24, 30, 29, 25,
                 20, 21, 32, 10, 31, 22, 33, 11, 32, 23, 34, 12, 33, 24, 35, 13, 34, 25, 18, 14};

  (*vert_adj_blocks) = {0,  2,  5,  8,  11, 14,  17,  21,  25,  29,  33, 36, 40,
                        44, 48, 52, 54, 57, 60,  63,  66,  68,  71,  75, 79, 83,
                        87, 89, 92, 95, 98, 101, 104, 108, 112, 116, 120};

  mesh_graph.build_from_offsets(std::move(vert_ids), std::move(vert_adj_blocks));

  graph::GraphPartitioner partitioner;
  std::vector<Uint> partition_id;

  partitioner.part_graph(mesh_graph, 4, "scotch", partition_id);

  for (auto id : partition_id)
  {
    std::cout << id << " ";
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
