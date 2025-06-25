#ifndef PDEKIT_Graph_Graph_Reordering_hpp
#define PDEKIT_Graph_Graph_Reordering_hpp

#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "graph/Graph.hpp"

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>

namespace pdekit
{

namespace graph
{

class GraphReordering
{
  public:
  // Default constructor
  GraphReordering();

  // Destructor
  ~GraphReordering();

  // Apply reverse Cuthill-McKee reordering through boost::graph library
  // algorithms
  template <typename VertType>
  static void compute_reverse_cuthill_mckee(const Graph<VertType> &graph,
                                            std::vector<VertType> &vertex_ids);

  // Apply reverse Cuthill-McKee reordering through boost::graph library
  // algorithms The input is a compressed graph
  template <typename VertType>
  static void compute_reverse_cuthill_mckee(const common::BlockArray<VertType, Uint> &crs_graph,
                                            std::vector<VertType> &vertex_ids);

  private:
};

// ----------------------------------------------------------------------------

template <typename VertType>
void GraphReordering::compute_reverse_cuthill_mckee(const Graph<VertType> &graph,
                                                    std::vector<VertType> &vertex_ids)
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                boost::property<boost::vertex_color_t, boost::default_color_type,
                                                boost::property<boost::vertex_degree_t, int>>>
      BoostGraph;

  typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<BoostGraph>::vertices_size_type size_type;

  // Create boost graph object and fill it with edges
  BoostGraph boost_graph(graph.nb_vertices());
  for (Uint v = 0; v < graph.nb_vertices(); ++v)
  {
    const std::pair<typename Graph<VertType>::adj_vertex_const_iterator,
                    typename Graph<VertType>::adj_vertex_const_iterator>
        adj_vert_iterators = graph.adjacent_vertices(v);

    for (typename Graph<VertType>::adj_vertex_const_iterator vert_it = adj_vert_iterators.first;
         vert_it != adj_vert_iterators.second; ++vert_it)
    {
      if (static_cast<VertType>(v) <= *vert_it)
      {
        boost::add_edge(v, *vert_it, boost_graph);
      }
    }
  }

  std::vector<Vertex> inv_perm(boost::num_vertices(boost_graph));
  // reverse cuthill-mckee ordering
  boost::cuthill_mckee_ordering(boost_graph, inv_perm.rbegin(),
                                boost::get(boost::vertex_color, boost_graph),
                                boost::make_degree_map(boost_graph));

  /*
  boost::property_map<BoostGraph, boost::vertex_index_t>::type
    index_map = boost::get(boost::vertex_index, boost_graph);
  */

  // The result is 'inv_perm', which is the permutation from the new ordering
  // to the old ordering.
  // We want the opposite:

  vertex_ids.resize(inv_perm.size());

  for (size_type i = 0; i < inv_perm.size(); ++i)
  {
    vertex_ids[inv_perm[i]] = i;
  }

  /*
  for (size_type i = 0; i != inv_perm.size(); ++i)
  {
    vertex_ids[index_map[inv_perm[i]]] = i;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename VertType>
void GraphReordering::compute_reverse_cuthill_mckee(
    const common::BlockArray<VertType, Uint> &crs_graph, std::vector<VertType> &vertex_ids)
{
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                                boost::property<boost::vertex_color_t, boost::default_color_type,
                                                boost::property<boost::vertex_degree_t, int>>>
      BoostGraph;

  typedef boost::graph_traits<BoostGraph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<BoostGraph>::vertices_size_type size_type;

  // Create boost graph object and fill it with edges
  BoostGraph boost_graph(crs_graph.nb_blocks());
  for (Uint v = 0; v < crs_graph.nb_blocks(); ++v)
  {
    common::ArrayView<const VertType, _1D, Uint> verts_adj_to_v = crs_graph.const_block(v);

    for (Uint ineighb = 0; ineighb < verts_adj_to_v.size(); ++ineighb)
    {
      if (static_cast<VertType>(v) <= verts_adj_to_v[ineighb])
      {
        boost::add_edge(v, verts_adj_to_v[ineighb], boost_graph);
      }
    }
  }

  std::vector<Vertex> inv_perm(boost::num_vertices(boost_graph));
  // reverse cuthill-mckee ordering
  boost::cuthill_mckee_ordering(boost_graph, inv_perm.rbegin(),
                                boost::get(boost::vertex_color, boost_graph),
                                boost::make_degree_map(boost_graph));

  /*
  boost::property_map<BoostGraph, boost::vertex_index_t>::type
    index_map = boost::get(boost::vertex_index, boost_graph);
  */

  // The result is 'inv_perm', which is the permutation from the new ordering
  // to the old ordering.
  // We want the opposite:

  vertex_ids.resize(inv_perm.size());

  for (size_type i = 0; i < inv_perm.size(); ++i)
  {
    vertex_ids[inv_perm[i]] = i;
  }

  /*
  for (size_type i = 0; i != inv_perm.size(); ++i)
  {
    vertex_ids[index_map[inv_perm[i]]] = i;
  }
  */
}

// ----------------------------------------------------------------------------

} // namespace graph

} // namespace pdekit

#endif
