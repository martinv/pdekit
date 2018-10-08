/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE pdekit_graph_utest
#include <boost/test/unit_test.hpp>

/// PDEKIT headers
#include "graph/GraphReordering.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(graph_constructor_utest)
{
  const Uint nb_vertices = 5;

  graph::Graph<Int> my_graph(nb_vertices);

  // Graph
  BOOST_CHECK_EQUAL(my_graph.nb_vertices(), nb_vertices);
  BOOST_CHECK_EQUAL(my_graph.nb_edges(), 0u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(graph_edge_insertion_utest)
{
  const Uint nb_vertices = 5;

  graph::Graph<Int> my_graph(nb_vertices);

  my_graph.insert_edge(2, 2);
  my_graph.insert_edge(2, 3);
  my_graph.insert_edge_unique(2, 3);
  my_graph.insert_edge(4, 1);

  // Graph
  BOOST_CHECK_EQUAL(my_graph.nb_vertices(), nb_vertices);
  BOOST_CHECK_EQUAL(my_graph.nb_edges(), 3u);

  common::BlockArray<Int, Uint> crs_graph;
  my_graph.compress_to_crs(crs_graph);

  BOOST_CHECK_EQUAL(crs_graph.nb_blocks(), 5u);
  const common::ArrayView<const Int, _1D, Uint> adj_verts_to0 = crs_graph.const_block(0);
  const common::ArrayView<const Int, _1D, Uint> adj_verts_to1 = crs_graph.const_block(1);
  const common::ArrayView<const Int, _1D, Uint> adj_verts_to2 = crs_graph.const_block(2);
  const common::ArrayView<const Int, _1D, Uint> adj_verts_to3 = crs_graph.const_block(3);
  const common::ArrayView<const Int, _1D, Uint> adj_verts_to4 = crs_graph.const_block(4);

  BOOST_CHECK_EQUAL(adj_verts_to0.size(), 0u);
  BOOST_CHECK_EQUAL(adj_verts_to1.size(), 0u);
  BOOST_CHECK_EQUAL(adj_verts_to2.size(), 2u);
  BOOST_CHECK_EQUAL(adj_verts_to3.size(), 0u);
  BOOST_CHECK_EQUAL(adj_verts_to4.size(), 1u);

  BOOST_CHECK_EQUAL(adj_verts_to2[0], 2u);
  BOOST_CHECK_EQUAL(adj_verts_to2[1], 3u);

  BOOST_CHECK_EQUAL(adj_verts_to4[0], 1u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(graph_adjacency_iterator_utest)
{
  const Uint nb_vertices = 5;

  graph::Graph<Int> my_graph(nb_vertices);
  typedef graph::Graph<Int>::adj_vertex_const_iterator adj_iterator_type;

  my_graph.insert_edge(2, 2);
  my_graph.insert_edge(2, 3);
  my_graph.insert_edge_unique(2, 3);
  my_graph.insert_edge(4, 1);
  my_graph.insert_edge(4, 3);
  my_graph.insert_edge(4, 2);

  std::pair<adj_iterator_type, adj_iterator_type> it_range = my_graph.adjacent_vertices(4);
  for (adj_iterator_type it = it_range.first; it != it_range.second; ++it)
  {
    std::cout << *it << " ";
  }
  std::cout << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(graph_reordering_utest)
{
  const Uint nb_vertices = 10;

  graph::Graph<Int> my_graph(nb_vertices);
  std::vector<graph::Graph<Int>::vertex_type> new_vertex_ids;

  my_graph.insert_edge_unique(0, 3);
  my_graph.insert_edge_unique(0, 5);
  my_graph.insert_edge_unique(1, 2);
  my_graph.insert_edge_unique(1, 4);
  my_graph.insert_edge_unique(1, 6);
  my_graph.insert_edge_unique(1, 9);
  my_graph.insert_edge_unique(2, 3);
  my_graph.insert_edge_unique(2, 4);
  my_graph.insert_edge_unique(3, 5);
  my_graph.insert_edge_unique(3, 8);
  my_graph.insert_edge_unique(4, 6);
  my_graph.insert_edge_unique(5, 6);
  my_graph.insert_edge_unique(5, 7);
  my_graph.insert_edge_unique(6, 7);

  graph::GraphReordering::compute_reverse_cuthill_mckee(my_graph, new_vertex_ids);
  BOOST_CHECK_EQUAL(new_vertex_ids.size(), nb_vertices);

  const std::vector<Int> vertex_id_check = {0, 8, 7, 4, 6, 2, 5, 3, 1, 9};
  for (Uint i = 0; i < vertex_id_check.size(); ++i)
  {
    BOOST_CHECK_EQUAL(new_vertex_ids[i], vertex_id_check[i]);
  }

  // my_graph.apply_reordering(new_vertex_ids);
}

// ----------------------------------------------------------------------------
