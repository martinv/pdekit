#ifndef PDEKIT_Graph_Graph_hpp
#define PDEKIT_Graph_Graph_hpp

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <vector>

#include "common/BlockArray.hpp"
#include "common/PDEKit.hpp"

namespace pdekit
{

namespace graph
{

template <typename VertType>
class Graph
{
  public:
  using vertex_type               = VertType;
  using adj_vertex_const_iterator = typename std::vector<VertType>::const_iterator;

  /// Constructor
  Graph(const Uint nb_vertices);

  /// Delete copy constructor
  Graph(const Graph &other) = delete;

  /// Destructor
  ~Graph();

  /// Deleted assignment operator
  Graph &operator=(const Graph &rhs) = delete;

  /// Get the number of vertices
  /// @return nb_vertices ... number of vertices in the graph
  Uint nb_vertices() const;

  /// Get the number of edges
  /// @return nb_edges ... number of edges in the graph
  Uint nb_edges() const;

  /// Add an oriented edge to graph, don't check
  /// whether the edge already exists in graph or not
  void insert_edge(const VertType vert_start, const VertType vert_end);

  /// Add an oriented edge to graph, AND check
  /// whether the edge already exists in graph or not
  /// Edge only added if it does not exist yet
  void insert_edge_unique(const VertType vert_start, const VertType vert_end);

  /// Remove all edges from the graph
  void remove_all_edges();

  /// Get a pair of constant iterators representing the half-closed range
  /// [w1,..., wn) which delimits all vertices adjacent to vertex v
  std::pair<adj_vertex_const_iterator, adj_vertex_const_iterator> adjacent_vertices(
      const VertType v) const;

  /// Get the number of adjacent vertices to vertex v
  Uint number_adj_vertices(const VertType v) const;

  /// Apply reordering of graph vertices
  /// new_vert_labels[i] contains the new index of graph vertex i
  void apply_reordering(const std::vector<VertType> &new_vert_labels);

  /// Store the graph values in CRS format
  template <typename SizeType>
  void compress_to_crs(common::BlockArray<VertType, SizeType> &crs_array) const;

  /// Remove duplicate neighbours of given vertex
  void remove_duplicate_neighbours(const VertType v);

  /// Write the graph to a simple text file
  void write_to_file(const std::string &filename) const;

  /// Write the graph to an svg file
  void write_to_svg(const std::string &filename) const;

  /// Print graph edges to screen
  void print() const;

  private:
  /// Graph stored as adjacency list
  ///  The size of m_graph corresponds to the number of vertices in graph
  /// m_graph[i] is a vector of all vertices connected to vertex i by an edge
  std::vector<std::unique_ptr<std::vector<VertType>>> m_graph;
};

// ----------------------------------------------------------------------------

template <typename VertType>
Graph<VertType>::Graph(const Uint nb_vertices)
{
  m_graph.resize(nb_vertices);
  for (Uint i = 0; i < m_graph.size(); ++i)
  {
    m_graph[i] = std::unique_ptr<std::vector<VertType>>(new std::vector<VertType>());
  }
}

// ----------------------------------------------------------------------------

template <typename VertType>
Graph<VertType>::~Graph()
{
}

// ----------------------------------------------------------------------------

template <typename VertType>
Uint Graph<VertType>::nb_vertices() const
{
  return m_graph.size();
}

// ----------------------------------------------------------------------------

template <typename VertType>
Uint Graph<VertType>::nb_edges() const
{
  Uint nb_edges = 0;

  for (Uint v = 0; v < m_graph.size(); ++v)
  {
    nb_edges += m_graph[v]->size();
  }
  return nb_edges;
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::insert_edge(const VertType vert_start, const VertType vert_end)
{
  m_graph[vert_start]->push_back(vert_end);
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::insert_edge_unique(const VertType vert_start, const VertType vert_end)
{
  std::vector<VertType> &node_edges = *m_graph[vert_start];

  /*
  typename std::vector<VertType>::iterator it =
      std::find(node_edges.begin(), node_edges.end(), vert_end);
  */
  auto it = std::find(node_edges.begin(), node_edges.end(), vert_end);

  if (it == node_edges.end())
  {
    node_edges.push_back(vert_end);
  }
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::remove_all_edges()
{
  for (Uint i = 0; i < m_graph.size(); ++i)
  {
    m_graph[i]->resize(0);
  }
}

// ----------------------------------------------------------------------------

template <typename VertType>
std::pair<typename Graph<VertType>::adj_vertex_const_iterator,
          typename Graph<VertType>::adj_vertex_const_iterator>
Graph<VertType>::adjacent_vertices(const VertType v) const
{
  return std::make_pair(m_graph[v]->cbegin(), m_graph[v]->cend());
}

// ----------------------------------------------------------------------------

template <typename VertType>
Uint Graph<VertType>::number_adj_vertices(const VertType v) const
{
  return m_graph[v]->size();
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::apply_reordering(const std::vector<VertType> &new_vert_labels)
{
  for (Uint v = 0; v < m_graph.size(); ++v)
  {
    std::vector<VertType> &adj_vertices = *m_graph[v];
    for (Uint i = 0; i < adj_vertices.size(); ++i)
    {
      adj_vertices[i] = new_vert_labels[adj_vertices[i]];
    }
  }

  std::vector<std::unique_ptr<std::vector<VertType>>> tmp(m_graph.size());
  for (Uint v = 0; v < m_graph.size(); ++v)
  {
    tmp[new_vert_labels[v]] = std::move(m_graph[v]);
  }

  std::swap(m_graph, tmp);
}

// ----------------------------------------------------------------------------

template <typename VertType>
template <typename SizeType>
void Graph<VertType>::compress_to_crs(common::BlockArray<VertType, SizeType> &crs_array) const
{
  std::unique_ptr<std::vector<VertType>> graph_edges(new std::vector<VertType>());
  std::unique_ptr<std::vector<SizeType>> graph_edge_lengths(new std::vector<SizeType>());

  Uint n_edges = 0;
  for (Uint r = 0; r < m_graph.size(); ++r)
  {
    if (m_graph[r])
    {
      const std::vector<VertType> &row = *(m_graph[r]);
      n_edges += row.size();
    }
  }

  graph_edges->resize(n_edges);
  graph_edge_lengths->resize(m_graph.size());

  Uint edge_pos = 0;

  for (Uint r = 0; r < m_graph.size(); ++r)
  {
    if (m_graph[r])
    {
      const std::vector<VertType> &row = *(m_graph[r]);

      for (Uint c = 0; c < row.size(); ++c)
      {
        (*graph_edges)[edge_pos++] = row[c];
      }
      (*graph_edge_lengths)[r] = row.size();
    }
    else
    {
      (*graph_edge_lengths)[r] = 0;
    }
  }

  crs_array.build(std::move(graph_edges), std::move(graph_edge_lengths));
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::remove_duplicate_neighbours(const VertType v)
{
  if (m_graph[v])
  {
    std::vector<VertType> &adj_neighbours = (*m_graph[v]);
    std::sort(adj_neighbours.begin(), adj_neighbours.end());
    auto last = std::unique(adj_neighbours.begin(), adj_neighbours.end());
    adj_neighbours.erase(last, adj_neighbours.end());
  }
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::write_to_file(const std::string &filename) const
{
  std::ofstream out_file;
  out_file.open(filename.c_str());
  out_file << m_graph.size() << std::endl;
  out_file << nb_edges() << std::endl;

  for (Uint node_start = 0; node_start < m_graph.size(); ++node_start)
  {
    const std::vector<VertType> &node_out_edges = *m_graph[node_start];
    for (Uint i = 0; i < node_out_edges.size(); ++i)
    {
      out_file << node_start << " " << node_out_edges[i] << std::endl;
    }
  }
  out_file.close();
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::write_to_svg(const std::string &filename) const
{
  VertType max_id_in_graph = VertType();

  for (Uint node_start = 0; node_start < m_graph.size(); ++node_start)
  {
    const std::vector<VertType> &node_out_edges = *m_graph[node_start];
    for (Uint i = 0; i < node_out_edges.size(); ++i)
    {
      max_id_in_graph = std::max(max_id_in_graph, node_out_edges[i]);
    }
  }

  const Uint m     = m_graph.size();
  const VertType n = max_id_in_graph + 1;

  std::ofstream out_file;
  out_file.open(filename.c_str());

  const Real border = 0.5;

  out_file << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" "
              "viewBox=\"0 0 "
           << n + border << " " << m + border
           << " \">\n"
              "<style type=\"text/css\" >\n"
              "     <![CDATA[\n"
              "      rect.pixel {\n"
              "          fill:   #000000;\n"
              "      }\n"
              "    ]]>\n"
              "  </style>\n\n"
              "   <rect width=\""
           << n + border << "\" height=\"" << m + border
           << "\" fill=\"rgb(128, 128, 128)\"/>\n"
              "   <rect x=\""
           << 0.5 * border << "\" y=\"" << 0.5 * border << "\" width=\"" << n << "\" height=\"" << m
           << "\" fill=\"rgb(255, 255, 255)\"/>\n\n";

  for (Uint node_start = 0; node_start < m_graph.size(); ++node_start)
  {
    const std::vector<VertType> &node_out_edges = *m_graph[node_start];
    for (Uint i = 0; i < node_out_edges.size(); ++i)
    {
      // The 0.05 offset is here to center the actual nonzero entries are
      // squares of size [0.9 x 0.9] . This means that the right and
      // bottom edge of the plot have a gap of width 1.0 - 0.9 = 0.1.
      // Therefore we divide this gap by 2 and offset all squares
      // representing nonzero entries by 0.1/2 = 0.05
      out_file << R"(  <rect class="pixel" x=")" << node_out_edges[i] + 0.5 * border + 0.05
               << "\" y=\"" << node_start + 0.5 * border + 0.05
               << "\" width=\".9\" height=\".9\"/>\n";
    }
  }

  out_file << "</svg>" << std::endl;
  out_file.close();
}

// ----------------------------------------------------------------------------

template <typename VertType>
void Graph<VertType>::print() const
{
  for (Uint node_start = 0; node_start < m_graph.size(); ++node_start)
  {
    const std::vector<VertType> &node_out_edges = *m_graph[node_start];
    for (Uint i = 0; i < node_out_edges.size(); ++i)
    {
      std::cout << "{" << node_start << "," << node_out_edges[i] << "} ";
    }
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace graph

} // namespace pdekit

#endif
