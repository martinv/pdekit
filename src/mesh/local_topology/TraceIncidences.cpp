#include "mesh/local_topology/TraceIncidences.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// CLASS IncidencePair
// ----------------------------------------------------------------------------

IncidenceEntry::IncidenceEntry() : cell_idx(INVALID_CELL_ID), local_id(INVALID_LOC_ENTITY_ID)
{
}

// ----------------------------------------------------------------------------

IncidenceEntry::IncidenceEntry(const Uint par_cell, const SUint loc_id)
    : cell_idx(par_cell), local_id(loc_id)
{
}

// ----------------------------------------------------------------------------
// CLASS LocalIncidences
// ----------------------------------------------------------------------------

TraceIncidences::TraceIncidences() : m_incidences(), m_permutations(), m_idx(INVALID_CELL_ID)
{
}

// ----------------------------------------------------------------------------

TraceIncidences::TraceIncidences(
    const common::ArrayView<const IncidenceEntry, _1D, Uint> &incidences,
    const common::ArrayView<const EntityDofRealign, _1D, Uint> &permutations, const Uint idx)
{
  m_incidences.resize(incidences.size());
  m_permutations.resize(permutations.size());

  for (Uint i = 0; i < incidences.size(); ++i)
  {
    m_incidences[i]   = incidences[i];
    m_permutations[i] = permutations[i];
  }
  m_idx = idx;
}

// ----------------------------------------------------------------------------

TraceIncidences::TraceIncidences(const TraceIncidences &other_inc)
{
  // NOTE THAT THIS CLASS IS ONLY __PROXY__ TO UNDERLYING DATA - WE DON'T
  // COPY THE DATA POINTED TO, BUT ONLY POINTERS. In this particular case,
  // this is desired behaviour
  m_incidences   = other_inc.m_incidences;
  m_permutations = other_inc.m_permutations;
  m_idx          = other_inc.m_idx;
}

// ----------------------------------------------------------------------------

TraceIncidences &TraceIncidences::operator=(const TraceIncidences &other_inc)
{
  // NOTE THAT THIS CLASS IS ONLY __PROXY__ TO UNDERLYING DATA - WE DON'T
  // COPY THE DATA POINTED TO, BUT ONLY POINTERS. In this particular case,
  // this is desired behaviour
  m_incidences   = other_inc.m_incidences;
  m_permutations = other_inc.m_permutations;
  m_idx          = other_inc.m_idx;
  return *this;
}

// ----------------------------------------------------------------------------

TraceIncidences::~TraceIncidences()
{
}

// ----------------------------------------------------------------------------

Uint TraceIncidences::size() const
{
  return m_incidences.size();
}

// ----------------------------------------------------------------------------

bool TraceIncidences::has_linear_cell_idx(const Uint cell_idx) const
{
  for (auto inc : m_incidences)
  {
    if (inc.cell_idx == cell_idx)
    {
      return true;
    }
  }
  return false;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const TraceIncidences &iarray)
{
  for (Uint i = 0; i < iarray.size(); ++i)
  {
    os << "[" << iarray.m_incidences[i].cell_idx << "," << iarray.m_incidences[i].local_id << "] "
       << " (" << iarray.m_permutations[i].get().code() << ")";
    if ((i + 1) < iarray.size())
    {
      os << "\n";
    }
  }
  /*
  os << "[" << (*(iarray.m_first_incidence + iarray.m_size - 1)).cell_idx <<
  ","
     << (*(iarray.m_first_incidence + iarray.m_size - 1)).local_id << "]";
  */
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
