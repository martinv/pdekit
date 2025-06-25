#ifndef PDEKIT_Mesh_Containers_Mesh_Statistics_hpp
#define PDEKIT_Mesh_Containers_Mesh_Statistics_hpp

#include <iostream>
#include <vector>

namespace pdekit
{

namespace mesh
{

class MeshStatistics
{

  public:
  /// Default constructor
  MeshStatistics();

  /// Copy constructor
  MeshStatistics(const MeshStatistics &other_stat);

  /// Assignment operator
  MeshStatistics &operator=(const MeshStatistics &other_stat);

  /// Default destructor
  ~MeshStatistics();

  public:
  int tot_nb_cells;
  int tot_nb_facets;

  std::vector<int> tot_nb_cells_in_level;
};

/// Output stream operator - print the statistics
std::ostream &operator<<(std::ostream &os, MeshStatistics const &stat);

} // namespace mesh

} // namespace pdekit

#endif
