#ifndef PDEKIT_Mesh_Containers_Tria_Facets_hpp
#define PDEKIT_Mesh_Containers_Tria_Facets_hpp

#include <chrono>
#include <forward_list>
#include <iostream>
#include <tuple>
#include <vector>

#include "mesh/CellBuffer.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/EntityStatus.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/TopologyPredicates.hpp"
#include "mesh/iterators/TraceTopologyIterator.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/view/TraceTopologyView.hpp"

namespace pdekit
{

namespace mesh
{

namespace internal
{

/**
 *  This class stores all facets in adaptive mesh. The facets are not stored
 * with respect to mesh levels. Neither are they ordered in any way except for
 * the following: Each facet can be seen from two (or more, in case of edges in
 * 3D) incident cells as a pair of indexes (cell_left, local_id_left) or
 * (cell_right, local_id_right). These two pairs are stored sequentially right
 * after each other for each mesh facet. This enables fast lookup of cells'
 * neighbours.
 *
 *  When the cell is refined, one facet (say the one on the left) might be
 * incident to multiple facets on the other side, say 2. In that case, the left
 * facet and the two right facets incident to it must again be stored
 * sequentially right after each other in one block. THE INCIDENT ENTITIES
 * FORMING ONE FACET BLOCK MUST 1) either belong to the same level 2) or their
 * level index may differ by 1 only (i.e. the block holds incidence information
 *     between facet of a cell C and the facets of children of some neighbor of
 * C)
 *
 *
 */

// Forward declaration to enable output operator
template <typename MeshConfig>
class TriaFacets;

template <typename MeshConfig>
std::ostream &operator<<(std::ostream &os, const TriaFacets<MeshConfig> &facets);

template <typename MeshConfig>
class TriaFacets
{
  public:
  /// TYPEDEFS
  using iterator = TraceTopologyIterator<TraceTopologyView<TriaFacets<MeshConfig>, ViewIsNotConst>>;
  using const_iterator =
      TraceTopologyIterator<TraceTopologyView<TriaFacets<MeshConfig>, ViewIsConst>>;

  enum
  {
    GDIM = MeshConfig::GDIM
  };

  enum
  {
    TDIM = MeshConfig::TDIM
  };

  enum
  {
    FacetDIM = MeshConfig::TDIM - 1
  };

  /// Constructor
  TriaFacets();

  /// Copy constructor, should not be used
  TriaFacets(const TriaFacets &other) = delete;

  /// Destructor
  ~TriaFacets();

  /// Assignement operator is disallowed
  TriaFacets &operator=(const TriaFacets &other) = delete;

  /// Return incidences corresponding to one facet, active or not
  const TraceIncidences facet_data(const FlatIdx facet_idx) const;

  /// Return incidences corresponding to one active facet
  const TraceIncidences active_facet_data(const ActiveIdx facet_idx) const;

  /// Return incidence arrays corresponding to one facet, active or not
  const std::tuple<common::ArrayView<const IncidenceEntry, _1D, Uint>,
                   common::ArrayView<const EntityDofRealign, _1D, Uint>>
  facet_view(const FlatIdx facet_idx) const;

  /// Return incidence arrays corresponding to one active facet
  const std::tuple<common::ArrayView<const IncidenceEntry, _1D, Uint>,
                   common::ArrayView<const EntityDofRealign, _1D, Uint>>
  active_facet_view(const ActiveIdx facet_idx) const;

  /// Number of facets, where each 'facet' can contain multiple pairs
  /// (cell_id, local_id). Only active facets are considered here
  Uint nb_active_facets() const;

  /// Number of facets stored in "TriaFacets"
  Uint nb_all_facets() const;

  /// Return total number of facet entries
  Uint nb_stored_entries() const;

  /// Get the dimension of these facets
  Uint dim() const;

  /// Merge these facets with other mesh facets, leaving the
  /// other facets empty
  void merge(TriaFacets &other_mf);

  /// @return an iterator pointing to the beginning of the facets,
  /// const version
  const_iterator begin() const;

  /// @return const iterator regardless of whether TriaFacets
  /// is a constant object or not
  const_iterator cbegin() const;

  /// @return an iterator pointing 1 position after the end of the cell array
  const_iterator end() const;

  /// @return const iterator regardless of whether TriaFacets
  /// is a constant object or not
  const_iterator cend() const;

  /// Clone the data
  static void clone(const TriaFacets &facets_in, TriaFacets &facets_out);

  /// Set the dimension of facets
  void set_dim(const Uint topo_dim);

  /// Clear the contents of this container
  void clear();

  /// Rebuild facets from raw data arrays
  /// This destroys all data that was previously allocated
  /// @note: this also destroys the input arrays!
  void create_from_raw_data(std::unique_ptr<std::vector<IncidenceEntry>> &&incidences,
                            std::unique_ptr<std::vector<EntityDofRealign>> &&facet_permutations,
                            std::unique_ptr<std::vector<Uint>> &&facet_data_offsets);

  /// Rebuild facets from cells stored in a buffer
  /// This destroys all data that was previously allocated
  void create_from_cell_buffer(const CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> &cell_buffer);

  /// Insert a new facet or multiple facets
  /// @param incidences          ... vector of incidences (pair [cell idx, local_id] expressing
  ///                                facets or edges)
  /// @param facet_permutations  ... vector of permutations of facets so that
  ///                                they match when multiple pairs [cell idx, local_id]
  ///                                express the same facet
  /// @note All parameter vectors above MUST HAVE THE SAME LENGTH

  const TraceIncidences add_facet(std::vector<IncidenceEntry> const &incidences,
                                  std::vector<EntityDofRealign> const &facet_permutations);

  /// Insert a new facet or multiple facets
  /// @param incidences          ... vector of incidences (pair [cell idx, local_id] expressing
  ///                                facets or edges)
  /// @param facet_permutations  ... vector of permutations of facets so that
  ///                                they match when multiple pairs [cell idx, local_id]
  ///                                express the same facet
  ///
  const TraceIncidences add_facet(
      const common::ArrayView<const IncidenceEntry, _1D, Uint> &incidences,
      const common::ArrayView<const EntityDofRealign, _1D, Uint> &facet_permutations);

  /// Remove some facets
  /// @param facet_ids ... vector of indices of facets that should be removed
  void remove_facets(const std::vector<FlatIdx> &facet_ids);

  /// Update the status of each facet based on whether the facet refers to
  /// active cells or not
  void update_status(const std::vector<EntityStatus> &cell_status);

  /// Update the status of each facet based on whether the facet refers to
  /// active cells or not
  void set_status(const EntityStatus status);

  /// Update cell ids
  /// @note This method updates the ids of cells that the facets refer to
  ///       This is necessary when the numbering of cells changes.
  /// @param new_cell_id ... a vector that holds on position [i] the new
  /// number of the i-th cell
  ///        The length of this vector should be at least equal to the maximum
  ///        value of 'cell_id' stored in m_incidences
  void update_cell_ids(const std::vector<FlatIdx> &new_cell_id);

  /// Overloaded output operator <<
  friend std::ostream &operator<<<MeshConfig>(std::ostream &os,
                                              const TriaFacets<MeshConfig> &facets);

  private:
  enum class FacetType : unsigned short
  {
    PartitionInternal   = 0,
    OnPartitionBoundary = 1,
    OnDomainBoundary    = 2
  };

  /// Fill incidences corresponding to one facet
  void fill_facet_data(const FlatIdx facet_idx, TraceIncidences &incidences) const;

  /// Merge two stl vectors. This is just a helper function
  /// @param V1 ... first vector, to which we will add vector V2
  /// @param V2 ... second vector, which will be merged into vector V2
  /// @note: both vectors will be modified: V1 will increase in size and V2
  /// will be
  ///         left empty
  template <typename T>
  void merge_stl_vectors(std::vector<T> &V1, std::vector<T> &V2);

  /// The topological dimension of these facets. If m_facet_dim is 2D,
  /// these facets are 2D faces in 3-dimensional mesh.
  /// If m_facet_dim is 1D, then these facets are edges in 2D mesh or edges in
  /// 3D mesh
  Uint m_facet_dim;

  /// Vector of incidences - consists of pairs (cell_id, local_id). Multiple
  /// incidences can be associated with one mesh facet.
  std::vector<IncidenceEntry> m_incidences;

  /// Position of active facets - this vector holds
  /// indexes of all facets in TriaFacets that are currently
  /// marked as active
  std::vector<Uint> m_active_facet_pos;

  /// Permutations of each entry of each facet. Has to have the same length
  /// as the vector of incidence pairs 'm_incidences'
  std::vector<EntityDofRealign> m_facet_permutations;

  /// Status of each facet: active/not active ...
  /// The size of this vector should be equal to m_facet_data_offsets.size() -
  /// 1, because the status should be unique to each block
  std::vector<std::tuple<EntityStatus, FacetType>> m_facet_info;

  /// Vector of offsets delimiting groups of incidence pairs that belong
  /// to the same facet
  std::vector<Uint> m_facet_data_offsets;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
TriaFacets<MeshConfig>::TriaFacets() : m_facet_dim(_0D)
{
  m_incidences.resize(0);
  m_active_facet_pos.resize(0);
  m_facet_permutations.resize(0);
  m_facet_info.resize(0);
  m_facet_data_offsets.resize(1);
  m_facet_data_offsets[0] = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
TriaFacets<MeshConfig>::~TriaFacets()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const TraceIncidences TriaFacets<MeshConfig>::facet_data(const FlatIdx facet_idx) const
{
  const Uint pos_idx = facet_idx.id();
  const Uint size    = m_facet_data_offsets[pos_idx + 1] - m_facet_data_offsets[pos_idx];
  const common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      m_incidences.data() + m_facet_data_offsets[pos_idx], size);
  const common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      m_facet_permutations.data() + m_facet_data_offsets[pos_idx], size);

  TraceIncidences block(incidences, permutations, pos_idx);
  return block;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const TraceIncidences TriaFacets<MeshConfig>::active_facet_data(
    const ActiveIdx facet_idx) const
{
  const Uint pos_idx = m_active_facet_pos[facet_idx.id()];
  const Uint size    = m_facet_data_offsets[pos_idx + 1] - m_facet_data_offsets[pos_idx];

  common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      &m_incidences[m_facet_data_offsets[pos_idx]], size);
  common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      &m_facet_permutations[m_facet_data_offsets[pos_idx]], size);

  TraceIncidences block(incidences, permutations, facet_idx.id());
  return block;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const std::tuple<common::ArrayView<const IncidenceEntry, _1D, Uint>,
                        common::ArrayView<const EntityDofRealign, _1D, Uint>>
TriaFacets<MeshConfig>::facet_view(const FlatIdx facet_idx) const
{
  const Uint pos_idx = facet_idx.id();
  const Uint size    = m_facet_data_offsets[pos_idx + 1] - m_facet_data_offsets[pos_idx];
  const common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      m_incidences.data() + m_facet_data_offsets[pos_idx], size);
  const common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      m_facet_permutations.data() + m_facet_data_offsets[pos_idx], size);

  return std::make_tuple(incidences, permutations);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const std::tuple<common::ArrayView<const IncidenceEntry, _1D, Uint>,
                        common::ArrayView<const EntityDofRealign, _1D, Uint>>
TriaFacets<MeshConfig>::active_facet_view(const ActiveIdx facet_idx) const
{
  const Uint pos_idx = m_active_facet_pos[facet_idx.id()];
  const Uint size    = m_facet_data_offsets[pos_idx + 1] - m_facet_data_offsets[pos_idx];

  common::ArrayView<const IncidenceEntry, _1D, Uint> incidences(
      &m_incidences[m_facet_data_offsets[pos_idx]], size);
  common::ArrayView<const EntityDofRealign, _1D, Uint> permutations(
      &m_facet_permutations[m_facet_data_offsets[pos_idx]], size);

  return std::make_tuple(incidences, permutations);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline Uint TriaFacets<MeshConfig>::nb_active_facets() const
{
  return m_active_facet_pos.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaFacets<MeshConfig>::nb_all_facets() const
{
  /*
  if (m_facet_data_offsets.size() <= 1)
  {
    return 0u;
  }
  return (m_facet_data_offsets.size() - 1);
  */

  return m_facet_info.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaFacets<MeshConfig>::nb_stored_entries() const
{
  return m_incidences.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaFacets<MeshConfig>::dim() const
{
  return m_facet_dim;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::merge(TriaFacets &other_mf)
{
  if (m_facet_dim != other_mf.m_facet_dim)
  {
    std::cerr << "TriaFacets::merge: trying to merge facets of different "
                 "dimensions - aborting!"
              << std::endl;
  }

  // If the other facets are empty, there is nothing to do ...
  if (other_mf.nb_all_facets() == 0)
  {
    return;
  }

  merge_stl_vectors(m_incidences, other_mf.m_incidences);
  merge_stl_vectors(m_facet_permutations, other_mf.m_facet_permutations);
  merge_stl_vectors(m_facet_info, other_mf.m_facet_info);

  const Uint nb_active_facets = m_active_facet_pos.size() + other_mf.m_active_facet_pos.size();
  m_active_facet_pos.resize(nb_active_facets);

  // Index of active facet
  Uint af = 0;
  for (Uint f = 0; f < m_facet_info.size(); ++f)
  {
    if (std::get<0>(m_facet_info[f]) == EntityStatus::Active)
    {
      m_active_facet_pos[af++] = f;
    }
  }

  const Uint last_old_offset = m_facet_data_offsets.back();
  const Uint old_size        = m_facet_data_offsets.size();
  // We will skip the element '0' in the offsets
  // 'other_mf.m_facet_data_offsets', hence '-1'
  const Uint new_size = m_facet_data_offsets.size() + other_mf.m_facet_data_offsets.size() - 1;
  m_facet_data_offsets.resize(new_size);

  for (Uint i = 1; i < other_mf.m_facet_data_offsets.size(); ++i)
  {
    m_facet_data_offsets[old_size + i - 1] = last_old_offset + other_mf.m_facet_data_offsets[i];
  }

  other_mf.m_facet_data_offsets.resize(1);
  other_mf.m_facet_data_offsets[0] = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename TriaFacets<MeshConfig>::const_iterator TriaFacets<MeshConfig>::begin() const
{
  return const_iterator(*this, 0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename TriaFacets<MeshConfig>::const_iterator TriaFacets<MeshConfig>::cbegin() const
{
  return const_iterator(*this, 0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename TriaFacets<MeshConfig>::const_iterator TriaFacets<MeshConfig>::end() const
{
  return const_iterator(*this, m_facet_info.size());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename TriaFacets<MeshConfig>::const_iterator TriaFacets<MeshConfig>::cend() const
{
  return const_iterator(*this, m_facet_info.size());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::set_dim(const Uint topo_dim)
{
  m_facet_dim = topo_dim;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::clone(const TriaFacets &facets_in, TriaFacets &facets_out)
{
  facets_out.m_facet_dim = facets_in.m_facet_dim;

  facets_out.m_incidences.resize(facets_in.m_incidences.size());
  std::copy(facets_in.m_incidences.begin(), facets_in.m_incidences.end(),
            facets_out.m_incidences.begin());

  facets_out.m_active_facet_pos.resize(facets_in.m_active_facet_pos.size());
  std::copy(facets_in.m_active_facet_pos.begin(), facets_in.m_active_facet_pos.end(),
            facets_out.m_active_facet_pos.begin());

  facets_out.m_facet_permutations.resize(facets_in.m_facet_permutations.size());
  std::copy(facets_in.m_facet_permutations.begin(), facets_in.m_facet_permutations.end(),
            facets_out.m_facet_permutations.begin());

  facets_out.m_facet_info.resize(facets_in.m_facet_info.size());
  std::copy(facets_in.m_facet_info.begin(), facets_in.m_facet_info.end(),
            facets_out.m_facet_info.begin());

  facets_out.m_facet_data_offsets.resize(facets_in.m_facet_data_offsets.size());
  std::copy(facets_in.m_facet_data_offsets.begin(), facets_in.m_facet_data_offsets.end(),
            facets_out.m_facet_data_offsets.begin());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::clear()
{
  m_incidences.resize(0);
  m_active_facet_pos.resize(0);
  m_facet_permutations.resize(0);
  m_facet_info.resize(0);
  m_facet_data_offsets.resize(1);
  m_facet_data_offsets[0] = 0;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::create_from_raw_data(
    std::unique_ptr<std::vector<IncidenceEntry>> &&incidences,
    std::unique_ptr<std::vector<EntityDofRealign>> &&facet_permutations,
    std::unique_ptr<std::vector<Uint>> &&facet_data_offsets)
{
  if (facet_data_offsets->size() < 2)
  {
    return;
  }

  const Uint nb_facets = facet_data_offsets->size() - 1;

  m_incidences.clear();
  m_incidences.swap(*incidences);
  m_active_facet_pos.resize(nb_facets);
  for (Uint i = 0; i < m_active_facet_pos.size(); ++i)
  {
    m_active_facet_pos[i] = i;
  }
  m_facet_permutations.clear();
  m_facet_permutations.swap(*facet_permutations);
  m_facet_info.resize(nb_facets);
  // m_status.assign(nb_facets, EntityStatus::Active);
  for (std::tuple<EntityStatus, FacetType> &f_info : m_facet_info)
  {
    std::get<0>(f_info) = EntityStatus::Active;
  }

  m_facet_data_offsets.clear();
  m_facet_data_offsets.swap(*facet_data_offsets);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::create_from_cell_buffer(
    const CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> &cell_buffer)
{
  // Cell connectivity is the connectivity which holds cells with highest
  // topological dimension in the mesh
  // In 2D, the cell connectivity has 2D cells (triangles, quads)
  // In 3D, cell connectivity has 3D cells (tetra, hexa, pyramids ...)
  // This is the input parameter 'master_cells'

  // Each edge/face (WHICH IS WHAT WE WANT TO COMPUTE HERE and this
  // corresponds to the input argument called 'subcells') is encoded by 2
  // numbers called 'IncidencePair': [cell_number, local_face_number].
  // Example: if the sub-cell is encoded as [10,2], it means it's the second
  // face/edge of the tenth (volume) cell. A volume cell is a triag/quad ...
  // in 2D and tetra/hexa/pyramid/prism ... in 3D
  //
  // The final container which we want to fill is called 'sub_cells' and it
  // will hold all newly computed entities - edges or faces
  //
  // For the computation, we also use a temporary hash vector called
  // 'tmp_subcells'. Each entry of 'tmp_subcells' holds
  // a list of lists of pairs [Uint,Uint]. In other words, 'tmp_subcells' is a
  //
  // vector < list <       list < IncidencePair >     >    >
  //                      [ THIS IS A MULTIEDGE ]
  //
  // (vector, where each entry is a list of lists and each 'inner' list holds
  // pairs of integers)
  //
  //  One list of pairs represents a 'multiedge/multiface': one face in the
  //  mesh
  // is shared by 2 cells, which means that it can be (for example)
  // expressed as [10,2] or [11,3] (the same face in the mesh is the second
  // local face of cell 10 or local face 3 of cell 11). If we computed edges
  // (instead of faces), the list can be longer than 2 - one edge in 3D mesh
  // can be shared by many elements. The vector 'tmp_subcells' defined below
  // has at each entry [i] a list of such 'multiedges/multifaces', because at
  // the position[i], we store all entities such that their smallest vertex
  // number is 'i'

  using TP = TopologyPredicates;

  const auto start = std::chrono::high_resolution_clock::now();

  std::cout << "TriaFacets " << m_facet_dim << "D: computing skeleton" << std::endl;

  // using subcell_list_type = std::list<std::pair<Uint,Uint> >;
  using IncidenceInfoTuple = std::tuple<Uint, Uint, Uint, EntityDofRealign>;
  using hash_vector_type   = std::unordered_map<Uint, std::vector<IncidenceInfoTuple>>;

  hash_vector_type tmp_subcells;

  // Variable to count how many skeleton edges/facets have been found
  // in the whole mesh
  // This variable does NOT count the multiplicity of multi-edges or
  // multi-faces
  Uint nb_skeleton_blocks = 0;

  // Count how many sub-cells (pairs [cell_id, local_id]) have been inserted
  // including their multiplicity
  // Example: if edge is shared by cells 104 and 207 (it is, say, the first
  // edge in cell 104 and fourth edge in cell 207), then it will be stored as
  // {[104,1],[207,4]} and counts as one inserted sub-entity but two inserted
  // pairs (entries)
  Uint nb_inserted_entries = 0;

  bool found_matching_subcell = false;

  // Loop over all cell groups in buffer

  EntityDofRealign p;

  // Loop over all cells
  for (Uint c = 0; c < cell_buffer.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cell_buffer.active_cell(ActiveIdx(c));

    for (Uint isub = 0; isub < cell.nb_sub_elements(m_facet_dim); ++isub)
    {
      MeshEntity sub_cell = cell.sub_entity(m_facet_dim, isub);

      const Uint sub_cell_hash = TP::min_entity_vertex(sub_cell);
      // const Uint sub_cell_hash = MeshEntityVertHasher::hash(sub_cell);
      // const Uint sub_cell_hash = TP::sum_min_max_entity_dof(sub_cell);

      // Loop over a list of entities which could possibly match the
      // entity 'sub_cell'
      std::vector<IncidenceInfoTuple> &subcell_array = tmp_subcells[sub_cell_hash];

      found_matching_subcell            = false;
      Uint idx_matching_entity_in_array = 0;

      // Make sure that NO operation resizing subcell_array such as
      // push_back(), emplace_back() etc. occurs INSIDE this for loop
      // It leads to bugs.
      for (Uint i = 0; i < subcell_array.size(); ++i)
      {
        const IncidenceInfoTuple id_tuple_in_list = subcell_array[i];

        // Get the permutation and its code
        const EntityRealignCode pcode = std::get<3>(id_tuple_in_list).get().code();

        // We consider ONLY entities whose permutation code is identity.
        // The methods TP::entities_match() and
        // TP::entities_match_reverse() are not designed to work for
        // cases where the first argument (first entity to match) is not
        // identity.
        //
        // If the subcell_array has nonzero size, there surely must
        // exist at least one such entity, because entities that didn't
        // find a matching counterpart, are added with identical
        // permutation
        //
        // This also saves work, as we skip checking unsuitable
        // candidates.

        if (pcode.nb_rotations() == 0 && pcode.nb_flips() == 0)
        {
          MeshEntity test_matching_sub_cell =
              cell_buffer.active_cell(ActiveIdx(std::get<0>(id_tuple_in_list)));
          test_matching_sub_cell.local_transform(m_facet_dim, std::get<1>(id_tuple_in_list));

          found_matching_subcell = TP::entities_match(test_matching_sub_cell, sub_cell, p);
          if (!found_matching_subcell)
          {
            found_matching_subcell =
                TP::entities_match_reverse(test_matching_sub_cell, sub_cell, p);
          }

          // If we found matching entity, no need to continue search:
          // break the loop
          if (found_matching_subcell)
          {
            idx_matching_entity_in_array = i;
            break;
          }

        } // If the (possibly) matching sub-cell has identical
          // permutation

      } // Loop over the list of (possibly) matching sub-cells

      if (found_matching_subcell)
      {
        subcell_array.emplace_back(cell.idx(), isub,
                                   std::get<2>(subcell_array[idx_matching_entity_in_array]), p);
        nb_inserted_entries++;
      }

      else //(found_matching_subcell)
      {
        EntityDofRealign p;
        p.change_type(sub_cell.pt_set_id(),
                      EntityRealignCode::identity(sub_cell.pt_set_id().elem_shape()));
        subcell_array.emplace_back(cell.idx(), isub, nb_skeleton_blocks, p);

        nb_skeleton_blocks++;
        nb_inserted_entries++;
      }

    } // Loop over all sub-entities of 'cell'

  } // iterator loop over all master cells

  std::cout << "Nb of skeleton blocks = " << nb_skeleton_blocks << std::endl;
  std::cout << "Nb of inserted entries = " << nb_inserted_entries << std::endl;

  // We store sub-cells including their multiplicity
  m_incidences.resize(nb_inserted_entries);

  m_active_facet_pos.resize(nb_skeleton_blocks);
  for (Uint i = 0; i < nb_skeleton_blocks; ++i)
  {
    m_active_facet_pos[i] = i;
  }

  m_facet_permutations.resize(nb_inserted_entries);
  m_facet_info.resize(nb_skeleton_blocks);
  // m_status.assign(m_status.size(), EntityStatus::Active);

  for (std::tuple<EntityStatus, FacetType> &f_info : m_facet_info)
  {
    std::get<0>(f_info) = EntityStatus::Active;
  }

  // Resize fact_block_offsets to size (nb_skeleton_items+1) and
  // assign 0 to all of its entries
  m_facet_data_offsets.assign(nb_skeleton_blocks + 1, 0u);

  // nb_skeleton_items = 0;
  // nb_inserted_entries = 0;

  // First calculate how many entities are in each block
  // and store this information in offsets
  for (auto it : tmp_subcells)
  {
    const std::vector<IncidenceInfoTuple> &subcell_array = it.second;
    for (Uint isub = 0; isub < subcell_array.size(); ++isub)
    {
      const Uint incidence_block_id = std::get<2>(subcell_array[isub]);
      m_facet_data_offsets[incidence_block_id + 1]++;
    }
  }

  // Now turn the data in offsets into real offsets
  for (Uint i = 0; i < (m_facet_data_offsets.size() - 1); ++i)
  {
    m_facet_data_offsets[i + 1] += m_facet_data_offsets[i];
  }

  // Helper array that indicates how many entries in given facet block
  // have been filled so far
  std::vector<Uint> next_free_idx_in_block(m_facet_data_offsets.size() - 1, 0);

  for (auto it : tmp_subcells)
  {
    // Loop over all incidence information and distribute it into correct
    // location in member arrays
    const std::vector<IncidenceInfoTuple> &subcell_array = it.second;
    for (Uint isub = 0; isub < subcell_array.size(); ++isub)
    {
      const Uint incidence_block_id = std::get<2>(subcell_array[isub]);
      const Uint entry_pos =
          m_facet_data_offsets[incidence_block_id] + next_free_idx_in_block[incidence_block_id];

      m_incidences[entry_pos].cell_idx = std::get<0>(subcell_array[isub]);
      m_incidences[entry_pos].local_id = std::get<1>(subcell_array[isub]);
      m_facet_permutations[entry_pos]  = std::get<3>(subcell_array[isub]);
      next_free_idx_in_block[incidence_block_id]++;
    }
  }

  const auto end                              = std::chrono::high_resolution_clock::now();
  const std::chrono::duration<double> elapsed = end - start;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "TriaFacets: building of mesh skeleton of dim " << m_facet_dim << " took "
            << elapsed.count() << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const TraceIncidences TriaFacets<MeshConfig>::add_facet(
    std::vector<IncidenceEntry> const &incidences,
    std::vector<EntityDofRealign> const &facet_permutations)
{
#if 0
  for (Uint i = 0; i < incidences.size(); ++i)
  {
    m_incidences.push_back(incidences[i]);
    m_facet_permutations.push_back(facet_permutations[i]);
  }

  // The facet is not listed as active (yet) in m_active_facet_pos, since its
  // status is 'pending'
  m_facet_info.push_back(
      std::make_tuple(EntityStatus::PendingRefinement, FacetType::PartitionInternal));

  const Uint new_facet_block_offset_start = m_facet_data_offsets.back();
  m_facet_data_offsets.push_back(new_facet_block_offset_start + incidences.size());

  const Uint facet_idx = m_facet_data_offsets.size() > 1 ? m_facet_data_offsets.size() - 2 : 0;

  common::ArrayView<const IncidenceEntry, _1D, Uint> inc_block(
      &m_incidences[m_facet_data_offsets[facet_idx]],
      m_facet_data_offsets[facet_idx + 1] - m_facet_data_offsets[facet_idx]);
  common::ArrayView<const EntityDofRealign, _1D, Uint> perm_block(
      &m_facet_permutations[m_facet_data_offsets[facet_idx]],
      m_facet_data_offsets[facet_idx + 1] - m_facet_data_offsets[facet_idx]);

  TraceIncidences block(inc_block, perm_block, facet_idx);

  return block;
#endif

  const common::ArrayView<const IncidenceEntry, _1D, Uint> incidences_view(incidences.data(),
                                                                           incidences.size());
  const common::ArrayView<const EntityDofRealign, _1D, Uint> facet_permutations_view(
      facet_permutations.data(), facet_permutations.size());
  return this->add_facet(incidences_view, facet_permutations_view);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const TraceIncidences TriaFacets<MeshConfig>::add_facet(
    const common::ArrayView<const IncidenceEntry, _1D, Uint> &incidences,
    const common::ArrayView<const EntityDofRealign, _1D, Uint> &facet_permutations)
{
  for (Uint i = 0; i < incidences.size(); ++i)
  {
    m_incidences.push_back(incidences[i]);
    m_facet_permutations.push_back(facet_permutations[i]);
  }

  // The facet is not listed as active (yet) in m_active_facet_pos, since its
  // status is 'pending'
  m_facet_info.push_back(
      std::make_tuple(EntityStatus::PendingRefinement, FacetType::PartitionInternal));

  const Uint new_facet_block_offset_start = m_facet_data_offsets.back();
  m_facet_data_offsets.push_back(new_facet_block_offset_start + incidences.size());

  const Uint facet_idx = m_facet_data_offsets.size() > 1 ? m_facet_data_offsets.size() - 2 : 0;

  common::ArrayView<const IncidenceEntry, _1D, Uint> inc_block(
      &m_incidences[m_facet_data_offsets[facet_idx]],
      m_facet_data_offsets[facet_idx + 1] - m_facet_data_offsets[facet_idx]);
  common::ArrayView<const EntityDofRealign, _1D, Uint> perm_block(
      &m_facet_permutations[m_facet_data_offsets[facet_idx]],
      m_facet_data_offsets[facet_idx + 1] - m_facet_data_offsets[facet_idx]);

  TraceIncidences block(inc_block, perm_block, facet_idx);

  return block;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::remove_facets(const std::vector<FlatIdx> &facet_ids)
{
  const Uint nb_old_facets =
      (m_facet_data_offsets.size() < 1) ? 0 : (m_facet_data_offsets.size() - 1);

  if (nb_old_facets == 0)
  {
    std::cerr << "TriaFacets::remove_facets: error, no facets to remove" << std::endl;
    return;
  }

  std::vector<bool> remove_facet(nb_old_facets);
  remove_facet.assign(nb_old_facets, false);

  for (auto &&id : facet_ids)
  {
    remove_facet[id.id()] = true;
  }

  // Position where data can be copied
  Uint fill_pos = 0;
  // Position to first entry that should be preserved
  Uint valid_data_pos = 0;

  Uint nb_new_facets = 0;

  // I) Modify incidence and permutation data
  Uint nb_entries_after_removal = 0;
  for (Uint f = 0; f < nb_old_facets; ++f)
  {
    const Uint block_data_len = m_facet_data_offsets[f + 1] - m_facet_data_offsets[f];
    // If facet data should be preserved and there were
    // some data removed prior to processing the facet,
    // copy the data on the first available position
    if (!remove_facet[f])
    {
      nb_entries_after_removal += block_data_len;
      if (fill_pos != m_facet_data_offsets[f])
      {
        for (Uint i = 0; i < block_data_len; ++i)
        {
          m_incidences[fill_pos]         = m_incidences[valid_data_pos];
          m_facet_permutations[fill_pos] = m_facet_permutations[valid_data_pos];
          fill_pos++;
          valid_data_pos++;
        }
      }
      else
      {
        fill_pos += block_data_len;
        valid_data_pos += block_data_len;
      }

      nb_new_facets++;
    }
    else
    {
      valid_data_pos += block_data_len;
    }
  }

  m_incidences.resize(nb_entries_after_removal);
  m_facet_permutations.resize(nb_entries_after_removal);

  fill_pos = 0;

  // II) Modify all remaining information
  // Turn temporarily offsets into block sizes
  // Size of block f is stored on position [f+1]
  // m_facet_data_offsets[0] will remain = 0
  for (Uint f = nb_old_facets; f > 0; --f)
  {
    m_facet_data_offsets[f] = m_facet_data_offsets[f] - m_facet_data_offsets[f - 1];
  }

  m_facet_data_offsets[0] = 0;

  for (Uint f = 0; f < nb_old_facets; ++f)
  {
    if (!remove_facet[f])
    {
      if (fill_pos != f)
      {
        m_facet_info[fill_pos]             = m_facet_info[f];
        m_facet_data_offsets[fill_pos + 1] = m_facet_data_offsets[f + 1];
      }
      fill_pos++;
    }
  }

  // Convert block sizes back to offsets
  for (Uint f = 1; f <= nb_new_facets; ++f)
  {
    m_facet_data_offsets[f] += m_facet_data_offsets[f - 1];
  }

  m_facet_info.resize(nb_new_facets);
  m_facet_data_offsets.resize(nb_new_facets + 1);

  // Update the positions of active facets
  m_active_facet_pos.resize(0);

  for (Uint f = 0; f < m_facet_info.size(); ++f)
  {
    if (std::get<0>(m_facet_info[f]) == EntityStatus::Active)
    {
      m_active_facet_pos.push_back(f);
    }
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::update_status(const std::vector<EntityStatus> &cell_status)
{
  std::cout << "UPDATING STATUS OF " << cell_status.size() << " CELLS AND " << m_facet_info.size()
            << " FACETS" << std::endl;
  TraceIncidences one_facet;
  // m_status.assign(m_status.size(), EntityStatus::Active);
  for (std::tuple<EntityStatus, FacetType> &f_info : m_facet_info)
  {
    std::get<0>(f_info) = EntityStatus::Active;
  }

  const Uint nb_all_facets = m_facet_info.size();
  Uint nb_active_facets    = m_facet_info.size();

  for (Uint f = 0; f < nb_all_facets; ++f)
  {
    fill_facet_data(FlatIdx(f), one_facet);

    // bool facet_is_active = true;

    for (Uint ib = 0; ib < one_facet.size(); ++ib)
    {
      if (cell_status[one_facet.cell_id(ib)] == EntityStatus::NotActive)
      {
        // facet_is_active = false;

        // If the facet is not marked as inactive yet, do it and
        // decrease the number of active facets Note that if we simply
        // set the facet as inactive every time we encounter a cell
        // which is incident to the facet and which is not active, we
        // might deactivate some facets twice (through the left AND
        // right cells)

        if (std::get<0>(m_facet_info[f]) != EntityStatus::NotActive)
        {
          std::get<0>(m_facet_info[f]) = EntityStatus::NotActive;
          nb_active_facets--;
        }
      } // if
    }   // loop over incidences within one facet block

    /*
    std::cout << "*********************************" << std::endl;
    std::cout << one_facet;
    std::cout << "Should be active: " << facet_is_active << std::endl
              << std::endl;
    */
  }

  m_active_facet_pos.resize(nb_active_facets);
  // Index of active facet
  Uint af = 0;
  for (Uint f = 0; f < m_facet_info.size(); ++f)
  {
    if (std::get<0>(m_facet_info[f]) == EntityStatus::Active)
    {
      m_active_facet_pos[af++] = f;
    }
  }

  /*
  std::cout << "***********************" << std::endl;
  std::cout << "Active facet numbering:" << std::endl;
  std::cout << "***********************" << std::endl;

  for (Uint f = 0; f < m_active_facet_pos.size(); ++f)
  {
    std::cout << "[Absolute(" << m_active_facet_pos[f] << "),Active(" << f <<
  ")]" << std::endl;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::set_status(const EntityStatus status)
{
  for (auto &facet_info : m_facet_info)
  {
    std::get<0>(facet_info) = status;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::update_cell_ids(const std::vector<FlatIdx> &new_cell_id)
{
  for (Uint i = 0; i < m_incidences.size(); ++i)
  {
    const Uint new_id        = new_cell_id[m_incidences[i].cell_idx].id();
    m_incidences[i].cell_idx = new_id;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaFacets<MeshConfig>::fill_facet_data(const FlatIdx facet_idx,
                                             TraceIncidences &incidences) const
{
  const Uint pos_idx = facet_idx.id();
  const Uint size    = m_facet_data_offsets[pos_idx + 1] - m_facet_data_offsets[pos_idx];
  common::ArrayView<const IncidenceEntry, _1D, Uint> inc_block(
      &m_incidences[m_facet_data_offsets[pos_idx]], size);
  common::ArrayView<const EntityDofRealign, _1D, Uint> perm_block(
      &m_facet_permutations[m_facet_data_offsets[pos_idx]], size);

  incidences = TraceIncidences(inc_block, perm_block, pos_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename T>
void TriaFacets<MeshConfig>::merge_stl_vectors(std::vector<T> &V1, std::vector<T> &V2)
{
  const Uint old_size = V1.size();
  const Uint new_size = V1.size() + V2.size();
  V1.resize(new_size);

  for (Uint i = 0; i < V2.size(); ++i)
  {
    V1[old_size + i] = V2[i];
  }

  V2.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::ostream &operator<<(std::ostream &os, const TriaFacets<MeshConfig> &facets)
{
  if (facets.m_facet_data_offsets.size() > 1)
  {
    for (Uint f = 0; f < (facets.m_facet_data_offsets.size() - 1); ++f)
    {
      os << "Facet block " << f << "       [" << std::get<0>(facets.m_facet_info[f]) << "]"
         << std::endl;
      const Uint block_begin = facets.m_facet_data_offsets[f];
      const Uint block_end   = facets.m_facet_data_offsets[f + 1];

      for (Uint ib = block_begin; ib < block_end; ++ib)
      {
        os << "  [" << facets.m_incidences[ib].cell_idx << "," << facets.m_incidences[ib].local_id
           << "] (" << facets.m_facet_permutations[ib].get().code().as_string() << ")" << std::endl;
      } // Loop over incidence entries in one block

    } // Loop over all facets
  }

  /*
  os << "***********************" << std::endl;
  os << "Active facet numbering:" << std::endl;
  os << "***********************" << std::endl;

  for (Uint f = 0; f < facets.m_active_facet_pos.size(); ++f)
  {
    os << "[Absolute(" << facets.m_active_facet_pos[f] << "),Active(" << f <<
  ")]" << std::endl;
  }
  */

  return os;
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
