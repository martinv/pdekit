#ifndef PDEKIT_Mesh_Trace_Incidences_hpp
#define PDEKIT_Mesh_Trace_Incidences_hpp

#include "common/ArrayView.hpp"
#include "mesh/EntityDofRealign.hpp"
#include "mesh/MeshConstants.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

struct IncidenceEntry
{
  /// Default constructor
  IncidenceEntry();

  /// Construct incidence from parent cell and local id
  IncidenceEntry(const Uint par_cell, const SUint loc_id);

  inline bool operator<(const IncidenceEntry &inc_r) const
  {
    if ((cell_idx == inc_r.cell_idx) && (local_id <= inc_r.local_id))
    {
      return true;
    }
    if (cell_idx < inc_r.cell_idx)
    {
      return true;
    }
    return false;
  }

  /// Id of the incident parent cell
  Uint cell_idx;

  /// Local id in the incident cell
  SUint local_id;
};

// ----------------------------------------------------------------------------

class TraceIncidences
{
  public:
  /// Default constructor
  TraceIncidences();

  /// Constructor: set the first cell incidence, first entity permutation,
  /// length and index of this incidence array (block)
  TraceIncidences(const common::ArrayView<const IncidenceEntry, _1D, Uint> &incidences,
                  const common::ArrayView<const EntityDofRealign, _1D, Uint> &permutations,
                  const Uint idx);

  /// Copy constructor
  TraceIncidences(const TraceIncidences &other_inc);

  /// Assignment operator
  TraceIncidences &operator=(const TraceIncidences &other_inc);

  /// Destructor
  ~TraceIncidences();

  /// Vector of all incidences
  const common::ArrayView<const IncidenceEntry, _1D, Uint> incidences() const;

  /// Vector of permutations
  const common::ArrayView<const EntityDofRealign, _1D, Uint> permutations() const;

  /// Return the index (number) of this incidence array
  Uint idx() const;

  /// Return the cell id of the i-th incidence in this object
  Uint cell_id(const Uint i) const;

  /// Return the local id of the i-th incidence
  SUint local_id(const Uint i) const;

  /// Return the cell id and local id as a pair of values
  IncidenceEntry incidence_entry(const Uint i) const;

  /// Return the permutation of the i-th incidence
  EntityDofRealign permutation(const Uint i) const;

  /// Return the number of entries stored
  Uint size() const;

  /// Return true if this incidence block contains given cell id
  bool has_linear_cell_idx(const Uint cell_idx) const;

  /// Output the contents of this incidence array
  friend std::ostream &operator<<(std::ostream &os, const TraceIncidences &iarray);

  private:
  /// Pointer to the first cell id (first incidence)
  std::vector<IncidenceEntry> m_incidences;

  /// Pointer to the permutation of the first cell
  std::vector<EntityDofRealign> m_permutations;

  /// Index (number) of this incidence.
  Uint m_idx;
};

// ----------------------------------------------------------------------------

inline const common::ArrayView<const IncidenceEntry, _1D, Uint> TraceIncidences::incidences() const
{
  return common::ArrayView<const IncidenceEntry, _1D, Uint>(m_incidences.data(),
                                                            m_incidences.size());
}

// ----------------------------------------------------------------------------

inline const common::ArrayView<const EntityDofRealign, _1D, Uint> TraceIncidences::permutations()
    const
{
  return common::ArrayView<const EntityDofRealign, _1D, Uint>(m_permutations.data(),
                                                              m_permutations.size());
}

// ----------------------------------------------------------------------------

inline Uint TraceIncidences::idx() const
{
  return m_idx;
}

// ----------------------------------------------------------------------------

inline Uint TraceIncidences::cell_id(const Uint i) const
{
  return m_incidences[i].cell_idx;
}

// ----------------------------------------------------------------------------

inline SUint TraceIncidences::local_id(const Uint i) const
{
  return m_incidences[i].local_id;
}

// ----------------------------------------------------------------------------

inline IncidenceEntry TraceIncidences::incidence_entry(const Uint i) const
{
  const IncidenceEntry ipair(m_incidences[i]);
  return ipair;
}

// ----------------------------------------------------------------------------

inline EntityDofRealign TraceIncidences::permutation(const Uint i) const
{
  return m_permutations[i];
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
