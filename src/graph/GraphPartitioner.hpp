#ifndef PDEKIT_Graph_Graph_Partitioner_hpp
#define PDEKIT_Graph_Graph_Partitioner_hpp

#include "PDEKit_Config.hpp"
#include "common/BlockArray.hpp"

#if PDEKIT_HAVE_SCOTCH
#include "scotch.h"
#endif

namespace pdekit
{

namespace graph
{

// ----------------------------------------------------------------------------

class GraphPartitioner
{
  public:
  /// Default constructor
  GraphPartitioner();

  /// Copy constructor
  GraphPartitioner(const GraphPartitioner &rhs);

  /// Default destructor
  ~GraphPartitioner();

  /// Assignment operator
  GraphPartitioner &operator=(const GraphPartitioner &rhs);

  /// Partition the graph
  /// @param graph        ... graph stored as block array. The indices in n-th
  /// block
  ///                         are all vertices which the n-th vertex is
  ///                         connected to
  /// @param nb_parts     ... number of parts in which the graph should be
  /// partitioned
  /// @param algorithm    ... type of algorithm to use for partitioning
  /// @param partition_id ... vector of indices denoting to which partition
  /// each vertex belongs
  template <typename VertType>
  void part_graph(const common::BlockArray<VertType, VertType> &graph, const Uint nb_parts,
                  const std::string &algorithm, std::vector<Uint> &partition_id) const;

  private:
};

// ----------------------------------------------------------------------------

template <typename VertType>
void GraphPartitioner::part_graph(const common::BlockArray<VertType, VertType> &graph,
                                  const Uint nb_parts, const std::string &algorithm,
                                  std::vector<Uint> &partition_id) const
{
#if PDEKIT_HAVE_SCOTCH
  std::cout << "Graph partitioner: going to partition graph with serial Scotch "
               "library."
            << std::endl;

  SCOTCH_Graph GraphSCOTCH;
  SCOTCH_Strat StratSCOTCH;

  SCOTCH_graphInit(&GraphSCOTCH);
  SCOTCH_stratInit(&StratSCOTCH);
  SCOTCH_Num baseval = 0;

  std::vector<SCOTCH_Num> edgetab; // Adjacency array which stores global indices
  std::vector<SCOTCH_Num> verttab; // Adjacency index array of size (vertnbr+1) - offsets

  edgetab.resize(graph.size());
  verttab.resize(graph.nb_blocks() + 1);
  verttab[0] = 0;

  Uint insert_idx = 0;

  for (Uint b = 0; b < graph.nb_blocks(); ++b)
  {
    const common::ArrayView<const VertType, _1D, VertType> one_vert_adj = graph.const_block(b);
    for (Uint v = 0; v < one_vert_adj.size(); ++v)
    {
      edgetab[insert_idx++] = static_cast<SCOTCH_Num>(one_vert_adj[v]);
    }
    verttab[b + 1] = verttab[b] + one_vert_adj.size();
  }

  // Adjacency end index array. By default, vendtab = verttab + 1
  SCOTCH_Num *vendtab = verttab.data() + 1;

  // Vertex label array
  SCOTCH_Num *vlbltab = nullptr;

  SCOTCH_Num edgenbr = edgetab.size();
  SCOTCH_Num vertnbr = verttab.size() - 1;

  // Arc load array. Must have size edgenbr if it exists
  SCOTCH_Num *edlotab = nullptr;

  // Vertex load array. Must have size edgenbr if it exists
  SCOTCH_Num *velotab = nullptr;

  // std::cout << "Number of edges in dual graph = " << edgenbr << std::endl;
  // std::cout << "Number of vertices in dual graph = " << vertnbr <<
  // std::endl;

  SCOTCH_graphBuild(&GraphSCOTCH, baseval, vertnbr, verttab.data(), vendtab, velotab, vlbltab,
                    edgenbr, edgetab.data(), edlotab);

  int ierr = SCOTCH_graphCheck(&GraphSCOTCH);

  if (ierr == 0)
  {
    std::cout << "SCOTCH Graph data is consistent" << std::endl;
  }
  else
  {
    std::cout << "SCOTCH Graph data is not consistent!" << std::endl;
  }

  // This vectors holds partition number for each vertex
  std::vector<SCOTCH_Num> part_id_tmp(vertnbr);

  SCOTCH_graphPart(&GraphSCOTCH, static_cast<SCOTCH_Num>(nb_parts), &StratSCOTCH,
                   part_id_tmp.data());

  SCOTCH_graphExit(&GraphSCOTCH);
  SCOTCH_stratExit(&StratSCOTCH);

  partition_id.resize(vertnbr);

  for (Uint i = 0; i < vertnbr; ++i)
  {
    partition_id[i] = static_cast<Uint>(part_id_tmp[i]);
  }

#else
  std::cerr << "Graph partitioner: can't partition graph because SCOTCH "
               "library was not found."
            << std::endl;
#endif
}

// ----------------------------------------------------------------------------

} // namespace graph

} // namespace pdekit

#endif
