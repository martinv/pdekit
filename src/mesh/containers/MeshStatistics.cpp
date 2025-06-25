#include <algorithm>

#include "common/PDEKit.hpp"
#include "mesh/containers/MeshStatistics.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

MeshStatistics::MeshStatistics() : tot_nb_cells(0), tot_nb_facets(0)
{
  tot_nb_cells_in_level.resize(0);
}

// ----------------------------------------------------------------------------

MeshStatistics::MeshStatistics(const MeshStatistics &other_stat)
{
  tot_nb_cells  = other_stat.tot_nb_cells;
  tot_nb_facets = other_stat.tot_nb_facets;
  tot_nb_cells_in_level.resize(other_stat.tot_nb_cells_in_level.size());
  std::copy(other_stat.tot_nb_cells_in_level.begin(), other_stat.tot_nb_cells_in_level.end(),
            tot_nb_cells_in_level.begin());
}

// ----------------------------------------------------------------------------

MeshStatistics &MeshStatistics::operator=(const MeshStatistics &other_stat)
{
  tot_nb_cells  = other_stat.tot_nb_cells;
  tot_nb_facets = other_stat.tot_nb_facets;
  std::copy(other_stat.tot_nb_cells_in_level.begin(), other_stat.tot_nb_cells_in_level.end(),
            tot_nb_cells_in_level.begin());
  return *this;
}

// ----------------------------------------------------------------------------

MeshStatistics::~MeshStatistics()
{
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, MeshStatistics const &stat)
{
  os << " ============ Mesh statistics ============ " << std::endl;
  os << " Total number of cells: " << stat.tot_nb_cells << std::endl;
  os << " Total number of cells per level: " << std::endl;
  for (Uint l = 0; l < stat.tot_nb_cells_in_level.size(); ++l)
  {
    os << "   Cells in level " << l << ": " << stat.tot_nb_cells_in_level[l] << std::endl;
  }

  os << " Total number of facets: " << stat.tot_nb_facets << std::endl;

  os << " ========================================= " << std::endl;
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
