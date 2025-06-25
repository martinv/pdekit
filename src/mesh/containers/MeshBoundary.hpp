#ifndef PDEKIT_Mesh_Mesh_Boundary_hpp
#define PDEKIT_Mesh_Mesh_Boundary_hpp

#include <queue>
#include <vector>

#include "common/IteratorRange.hpp"
#include "common/PtrHandle.hpp"
#include "mesh/MeshConfig.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/containers/TriaFacets.hpp"
#include "mesh/iterators/BdryDofIterator.hpp"
#include "mesh/iterators/TraceTopologyIterator.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/view/BdryDofView.hpp"
#include "mesh/view/CellTopologyView.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

template <typename MeshConfig>
class DofMap;

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
class BoundaryFacets
{
  public:
  using mesh_cells_t = internal::TriaCells<MeshConfig>;

  using dof_iterator =
      BdryDofIterator<BdryDofView<DofMap<MeshConfig>, BoundaryFacets, ViewIsNotConst>>;
  using const_dof_iterator =
      BdryDofIterator<BdryDofView<DofMap<MeshConfig>, BoundaryFacets, ViewIsConst>>;

  using dof_iterator_typed =
      BdryDofIterator<BdryDofView<DofMap<MeshConfig>, BoundaryFacets, ViewIsNotConst>,
                      DofIterFilterTyped>;
  using const_dof_iterator_typed =
      BdryDofIterator<BdryDofView<DofMap<MeshConfig>, BoundaryFacets, ViewIsConst>,
                      DofIterFilterTyped>;

  using const_bdry_topo_iterator =
      TraceTopologyIterator<TraceTopologyView<internal::TriaFacets<MeshConfig>, ViewIsConst>>;

  /*
  using const_dof_range_typed =
  common::IteratorRange<const_dof_iterator_typed>; using all_dof_type_ranges =
  std::vector<const_dof_range_typed>;
  */

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
    BcDIM = BcDim
  };

  /// Default constructor
  BoundaryFacets(const common::PtrHandle<mesh_cells_t const> mesh_cells, const std::string &name);

  /// Copy constructor is deleted
  BoundaryFacets(const BoundaryFacets &other) = delete;

  /// Default destructor
  ~BoundaryFacets() = default;

  /// Assignment operator is deleted
  BoundaryFacets &operator=(const BoundaryFacets &rhs) = delete;

  /// Get the dimension of these facets
  constexpr Uint dim() const;

  /// Get the name of these boundary facets
  const std::string &name() const;

  /// Return the material id of this boundary
  Uint material_id() const;

  /// Return the topology cell adjacent to the i-th facet
  const CellTopologyView<MeshConfig> adjacent_tcell(const ActiveIdx facet_id) const;

  /// Return the information about one boundary cell. This info is a pair
  /// [cell index, local id], which says to which cell does the facet belong
  /// (cell index) and what is the local number of the facet in this cell.
  /// @param: entry_idx number of the facet in this boundary facet group
  // const IncidencePair bdry_cell_id(const Uint entry_idx) const;

  /// Return the information about one boundary cell.
  /// @param: idx number of the facet in this boundary facet group
  const MeshEntity active_cell(const DofMap<MeshConfig> &dof_handler, const ActiveIdx idx) const;

  /// Return the information about volume cell which is touching
  /// a boundary cell
  /// @param: idx number of the facet in this boundary facet group
  const MeshEntity adjacent_active_vol_cell(const DofMap<MeshConfig> &dof_handler,
                                            const ActiveIdx idx) const;

  /// Return the index of (volume) cell this facet belongs to
  /// @param  facet_id ... active index of facet on boundary
  /// @return flat (linearized) index of volume cell that contains
  ///         that facet
  FlatIdx parent_cell_id(const ActiveIdx facet_id) const;

  /// Local index in parent (adjacent volume cell)
  /// @param: idx number of the facet in this boundary facet group
  Uint local_id(const ActiveIdx facet_id) const;

  /// Fill one active cell
  void fill_active_cell(const DofMap<MeshConfig> &dof_handler, const ActiveIdx idx,
                        MeshEntity &bdry_cell) const;

  /// Return the number of facets/edges
  Uint nb_active_cells() const;

  /// Fill the indexes of facets that lie on the boundary
  void emplace_bdry_cell_ids(std::unique_ptr<std::vector<IncidenceEntry>> &&facet_ids,
                             const Uint material_id);

  /// Clone all data EXCEPT the pointer to cells
  static void clone_data(BoundaryFacets const &facets_in, BoundaryFacets &facets_out);

  /// Update the facet indices after a topological change (refinement) in
  /// MeshTopology
  void update();

  /// Return const iterator regardless of whether BoundaryFacets is a constant
  /// object or not
  const_dof_iterator cbegin(DofMap<MeshConfig> const &dof_handler) const;

  /// Return const iterator regardless of whether BoundaryFacets is a constant
  /// object or not
  const_dof_iterator cend(DofMap<MeshConfig> const &dof_handler) const;

  /// Return const iterator to the beginning of this boundary
  const_bdry_topo_iterator cbegin() const;

  /// Return const iterator to the end of this boundary
  const_bdry_topo_iterator cend() const;

  /// Fill vector of ranges for each element type
  /// There will be as many ranges as there are element types on the boundary
  /// of given dof_handler
  /// @param dof_handler ... dof handler whose boundary should be inspected
  /// @param ranges ... (output) vector of ranges, each range corresponds to
  /// the
  ///                   iterator pair [first,last) for one element type on the
  ///                   boundary
  void all_bdry_dof_ranges(
      DofMap<MeshConfig> const &dof_handler,
      std::vector<common::IteratorRange<const_dof_iterator_typed>> &ranges) const;

  private:
  /// Pointer to cell data storage
  common::PtrHandle<mesh_cells_t const> m_mesh_cells;

  /// Name of this group of boundary facets
  std::string m_name;

  /// Physical tag/material id
  Uint m_material_id;

  /// Vector of connectivity data: for each boundary cell,
  /// the absolute position/linear index of that cell and its boundary
  /// face index are stored
  /// This vector should be always updated so that it contains
  /// data for ACTIVE CELLS ONLY
  internal::TriaFacets<MeshConfig> m_active_facets;

  /// Facets on level zero (when the boundary facets are created)
  internal::TriaFacets<MeshConfig> m_zero_level_facets;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
BoundaryFacets<MeshConfig, BcDim>::BoundaryFacets(
    const common::PtrHandle<const mesh_cells_t> mesh_cells, const std::string &name)
    : m_mesh_cells(mesh_cells), m_name(name), m_material_id(0u)
{
  m_active_facets.set_dim(BcDim);
  m_zero_level_facets.set_dim(BcDim);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
constexpr Uint BoundaryFacets<MeshConfig, BcDim>::dim() const
{
  return BcDim;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
const std::string &BoundaryFacets<MeshConfig, BcDim>::name() const
{
  return m_name;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
Uint BoundaryFacets<MeshConfig, BcDim>::material_id() const
{
  return m_material_id;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
const CellTopologyView<MeshConfig> BoundaryFacets<MeshConfig, BcDim>::adjacent_tcell(
    const ActiveIdx facet_id) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(facet_id);
  return (*m_mesh_cells).cell(FlatIdx(bdry_facet.cell_id(0)));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
inline const MeshEntity BoundaryFacets<MeshConfig, BcDim>::active_cell(
    const DofMap<MeshConfig> &dof_handler, const ActiveIdx idx) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(idx);
  const IncidenceEntry inc_entry   = bdry_facet.incidence_entry(0);
  MeshEntity bdry_cell             = dof_handler.cell(FlatIdx(inc_entry.cell_idx));
  bdry_cell.local_transform(MeshConfig::TDIM - 1, inc_entry.local_id);
  bdry_cell.set_idx(idx.id());
  return bdry_cell;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
inline const MeshEntity BoundaryFacets<MeshConfig, BcDim>::adjacent_active_vol_cell(
    const DofMap<MeshConfig> &dof_handler, const ActiveIdx idx) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(idx);
  MeshEntity bdry_cell             = dof_handler.cell(FlatIdx(bdry_facet.cell_id(0)));
  return bdry_cell;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
inline FlatIdx BoundaryFacets<MeshConfig, BcDim>::parent_cell_id(const ActiveIdx facet_id) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(facet_id);
  return FlatIdx(bdry_facet.cell_id(0));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
inline Uint BoundaryFacets<MeshConfig, BcDim>::local_id(const ActiveIdx facet_id) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(facet_id);
  return bdry_facet.local_id(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
inline void BoundaryFacets<MeshConfig, BcDim>::fill_active_cell(
    const DofMap<MeshConfig> &dof_handler, const ActiveIdx idx, MeshEntity &bdry_cell) const
{
  const TraceIncidences bdry_facet = m_active_facets.active_facet_data(idx);
  const IncidenceEntry inc_entry   = bdry_facet.incidence_entry(0);
  dof_handler.fill_cell(FlatIdx(inc_entry.cell_idx), bdry_cell);
  bdry_cell.local_transform(MeshConfig::TDIM - 1, inc_entry.local_id);
  bdry_cell.set_idx(idx.id());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
Uint BoundaryFacets<MeshConfig, BcDim>::nb_active_cells() const
{
  return m_active_facets.nb_active_facets();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
void BoundaryFacets<MeshConfig, BcDim>::emplace_bdry_cell_ids(
    std::unique_ptr<std::vector<IncidenceEntry>> &&facet_ids, const Uint material_id)
{
  const Uint nb_bdry_facets = facet_ids->size();
  std::unique_ptr<std::vector<EntityDofRealign>> facet_permutations(
      new std::vector<EntityDofRealign>());
  facet_permutations->resize(nb_bdry_facets);

  for (Uint bc = 0; bc < nb_bdry_facets; ++bc)
  {
    const IncidenceEntry ientry              = (*facet_ids)[bc];
    const CellTopologyView<MeshConfig> tcell = (*m_mesh_cells).cell(FlatIdx(ientry.cell_idx));
    const std::shared_ptr<const StdRegionEntity> facet =
        tcell.sub_entity(MeshConfig::TDIM - 1, ientry.local_id);

    const PointSetTag facet_tag = facet->pt_set_id();
    const ElemShape facet_shape = facet_tag.elem_shape();
    const EntityRealignCode facet_realign_code(facet_shape, CellTransform::NO_TRANS, 0, facet_shape,
                                               0, 0);
    (*facet_permutations)[bc].change_type(std::make_pair(facet_tag, facet_realign_code));
  }

  std::unique_ptr<std::vector<Uint>> facet_data_offsets(new std::vector<Uint>());
  facet_data_offsets->resize(nb_bdry_facets + 1);
  (*facet_data_offsets)[0] = 0;
  for (Uint i = 1; i < facet_data_offsets->size(); ++i)
  {
    (*facet_data_offsets)[i] = (*facet_data_offsets)[i - 1] + 1;
  }

  m_active_facets.create_from_raw_data(std::move(facet_ids), std::move(facet_permutations),
                                       std::move(facet_data_offsets));
  internal::TriaFacets<MeshConfig>::clone(m_active_facets, m_zero_level_facets);

  m_material_id = material_id;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
void BoundaryFacets<MeshConfig, BcDim>::clone_data(BoundaryFacets const &facets_in,
                                                   BoundaryFacets &facets_out)
{
  facets_out.m_name        = facets_in.m_name;
  facets_out.m_material_id = facets_in.m_material_id;

  internal::TriaFacets<MeshConfig>::clone(facets_in.m_active_facets, facets_out.m_active_facets);
  internal::TriaFacets<MeshConfig>::clone(facets_in.m_zero_level_facets,
                                          facets_out.m_zero_level_facets);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
void BoundaryFacets<MeshConfig, BcDim>::update()
{
#if 0
  std::vector<IncidencePair> new_id_linear_idx;
  new_id_linear_idx.reserve(m_facet_id_linear_idx.size());
  new_id_linear_idx.resize(0);

  for (Uint c = 0; c < m_facet_id_linear_idx.size(); ++c)
  {
    const CellTopologyView<MeshConfig> tcell = (*m_mesh).cell(m_facet_id_linear_idx[c].cell_idx);
    if (tcell.status() == EntityStatus::NotActive)
    {
      const adapt::CellAdaptOp adapt_op = tcell.cell_adapt_op();
      const LocalIncidences facet_incidences =
          adapt_op.get().parent_child_incidences_on_facet(m_facet_id_linear_idx[c].local_id);

      const std::vector<CellTopologyView<MeshConfig>> cell_children = tcell.children();

      for (Uint inc = 0; inc < facet_incidences.size(); ++inc)
      {
        const Uint child_cell_id = cell_children[facet_incidences.cell_id(inc)].linear_pos_idx();
        const Uint local_facet_id = facet_incidences.local_id(inc);
        new_id_linear_idx.push_back(IncidencePair(child_cell_id, local_facet_id));
      }
    }
    else
    {
      new_id_linear_idx.push_back(m_facet_id_linear_idx[c]);
    }
  }

  m_facet_id_linear_idx.swap(new_id_linear_idx);
#endif

  // --------------------------------------------------------------------------------------------

  std::queue<IncidenceEntry> facets_to_check;
  std::unique_ptr<std::vector<IncidenceEntry>> facet_id_linear_idx(
      new std::vector<IncidenceEntry>());

  // First register all active cells on refinement level 0
  // If a cell is on level 0 but not active, then its children might be active
  // In this case, insert the level in queue
  for (Uint i = 0; i < m_zero_level_facets.nb_active_facets(); ++i)
  {
    const TraceIncidences zero_level_incidences =
        m_zero_level_facets.active_facet_data(ActiveIdx(i));

    const CellTopologyView<MeshConfig> tcell =
        (*m_mesh_cells).cell(FlatIdx(zero_level_incidences.cell_id(0)));
    if (tcell.refinement_level() == 0)
    {
      if (tcell.status() == EntityStatus::Active)
      {
        facet_id_linear_idx->push_back(zero_level_incidences.incidence_entry(0));
      }
      else
      {
        facets_to_check.push(zero_level_incidences.incidence_entry(0));
      }
    }
  }

  while (!facets_to_check.empty())
  {
    const IncidenceEntry facet_id = facets_to_check.front();
    facets_to_check.pop();

    const CellTopologyView<MeshConfig> tcell = (*m_mesh_cells).cell(FlatIdx(facet_id.cell_idx));
    const adapt::CellAdaptOp adapt_op        = tcell.cell_adapt_op();

    const TraceIncidences facet_incidences =
        adapt_op.get().parent_child_incidences_on_facet(facet_id.local_id);

    const std::vector<CellTopologyView<MeshConfig>> cell_children = tcell.children();

    for (Uint inc = 0; inc < facet_incidences.size(); ++inc)
    {
      const CellTopologyView<MeshConfig> child_tcell = cell_children[facet_incidences.cell_id(inc)];
      const Uint child_cell_id                       = child_tcell.linear_pos_idx().id();
      const SUint local_facet_id                     = facet_incidences.local_id(inc);

      if (child_tcell.status() == EntityStatus::Active)
      {
        facet_id_linear_idx->push_back(IncidenceEntry(child_cell_id, local_facet_id));
      }
      else
      {
        facets_to_check.push(IncidenceEntry(child_cell_id, local_facet_id));
      }
    } // Loop over facet_incidences
  }   // while the queue is not empty

  const Uint nb_bdry_facets = facet_id_linear_idx->size();
  std::unique_ptr<std::vector<EntityDofRealign>> facet_permutations(
      new std::vector<EntityDofRealign>());
  facet_permutations->resize(nb_bdry_facets);

  for (Uint bc = 0; bc < nb_bdry_facets; ++bc)
  {
    const IncidenceEntry ientry              = (*facet_id_linear_idx)[bc];
    const CellTopologyView<MeshConfig> tcell = (*m_mesh_cells).cell(FlatIdx(ientry.cell_idx));
    const std::shared_ptr<const StdRegionEntity> facet =
        tcell.sub_entity(MeshConfig::TDIM - 1, ientry.local_id);

    const PointSetTag facet_tag = facet->pt_set_id();
    const ElemShape facet_shape = facet_tag.elem_shape();
    const EntityRealignCode facet_realign_code(facet_shape, CellTransform::NO_TRANS, 0, facet_shape,
                                               0, 0);
    (*facet_permutations)[bc].change_type(std::make_pair(facet_tag, facet_realign_code));
  }

  std::unique_ptr<std::vector<Uint>> facet_data_offsets(new std::vector<Uint>());
  facet_data_offsets->resize(nb_bdry_facets + 1);
  (*facet_data_offsets)[0] = 0;
  for (Uint i = 1; i < facet_data_offsets->size(); ++i)
  {
    (*facet_data_offsets)[i] = (*facet_data_offsets)[i - 1] + 1;
  }

  m_active_facets.create_from_raw_data(
      std::move(facet_id_linear_idx), std::move(facet_permutations), std::move(facet_data_offsets));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
typename BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator BoundaryFacets<
    MeshConfig, BcDim>::cbegin(DofMap<MeshConfig> const &dof_handler) const
{
  const_dof_iterator it(dof_handler, *this, ActiveIdx(0));
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
typename BoundaryFacets<MeshConfig, BcDim>::const_dof_iterator BoundaryFacets<
    MeshConfig, BcDim>::cend(DofMap<MeshConfig> const &dof_handler) const
{
  const_dof_iterator it(dof_handler, *this, ActiveIdx(m_active_facets.nb_active_facets()));
  return it;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
typename BoundaryFacets<MeshConfig, BcDim>::const_bdry_topo_iterator BoundaryFacets<
    MeshConfig, BcDim>::cbegin() const
{
  return const_bdry_topo_iterator(m_active_facets, 0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, Uint BcDim>
typename BoundaryFacets<MeshConfig, BcDim>::const_bdry_topo_iterator BoundaryFacets<
    MeshConfig, BcDim>::cend() const
{
  return const_bdry_topo_iterator(m_active_facets, m_active_facets.nb_active_facets());
}

// -------------------------------------------------------------------------`---

template <typename MeshConfig, Uint BcDim>
void BoundaryFacets<MeshConfig, BcDim>::all_bdry_dof_ranges(
    DofMap<MeshConfig> const &dof_handler,
    std::vector<common::IteratorRange<const_dof_iterator_typed>> &ranges) const
{
  // The map 'iter_range' stores the newly found ranges for each element type
  // on boundary
  using range_map =
      std::map<PointSetTag, std::pair<const_dof_iterator_typed, const_dof_iterator_typed>>;
  range_map iter_ranges;
  for (Uint c = 0; c < m_active_facets.nb_active_facets(); ++c)
  {
    const TraceIncidences bdry_facet = m_active_facets.active_facet_data(ActiveIdx(c));

    // Get the type of current boundary cell
    MeshEntity facet = dof_handler.cell(FlatIdx(bdry_facet.cell_id(0)));
    facet.local_transform(MeshConfig::TDIM - 1, bdry_facet.local_id(0));
    const PointSetTag cell_tag = facet.pt_set_id();
    // Check if range corresponding to the cell type is already present in
    // the map If not, insert a new range, else modify the 'end' iterator
    // Note that at this moment, the 'end' iterator points to the last
    // element of given type, NOT one element past!
    const typename range_map::iterator map_iter = iter_ranges.find(cell_tag);
    if (map_iter == iter_ranges.end())
    {
      // const_dof_iterator_typed new_iter_typed_begin(dof_handler, *this,
      // cell_tag, 0);
      const_dof_iterator_typed new_iter_typed_begin(dof_handler, *this, ActiveIdx(c), cell_tag);
      const_dof_iterator_typed new_iter_typed_end = new_iter_typed_begin;
      ++new_iter_typed_end;

      // synchronize_dof_iterators(all_dof_it, new_iter_typed_begin);
      std::pair<const_dof_iterator_typed, const_dof_iterator_typed> new_iter_range(
          new_iter_typed_begin, new_iter_typed_end);
      iter_ranges.insert(
          std::pair<PointSetTag, std::pair<const_dof_iterator_typed, const_dof_iterator_typed>>(
              cell_tag, new_iter_range));
    }
    else
    {
      // Else just modify the 'end' iterator in the range
      // const_dof_iterator_typed new_iter_typed(dof_handler, *this,
      // cell_tag, 0); synchronize_dof_iterators(all_dof_it,
      // new_iter_typed);

      const_dof_iterator_typed new_iter_typed(dof_handler, *this, ActiveIdx(c), cell_tag);
      ++new_iter_typed;
      map_iter->second.second = new_iter_typed;
    }
  }

  // Clear the original vector of ranges
  ranges.clear();
  for (typename range_map::const_iterator map_iter = iter_ranges.cbegin();
       map_iter != iter_ranges.cend(); ++map_iter)
  {
    ranges.push_back(common::IteratorRange<const_dof_iterator_typed>(map_iter->second.first,
                                                                     map_iter->second.second));
  }
}

// ----------------------------------------------------------------------------
// This class holds a group of domains representing the mesh boundaries
// ----------------------------------------------------------------------------

template <typename MeshConfig>
class MeshBoundarySet
{
  public:
  using mesh_cells_t = internal::TriaCells<MeshConfig>;

  using bdry_facets_shared_ptr = std::shared_ptr<BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>>;

  using const_bdry_topo_iterator =
      typename BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>::const_bdry_topo_iterator;

  /// Construct a SubCells object
  MeshBoundarySet();

  /// Destroy SubCellss
  ~MeshBoundarySet();

  void init(const mesh_cells_t &mesh_cells);

  /// Reset all data
  void clear();

  /// Copy the contents of one container into another
  /// Note that only the data in Cells is copied, but
  /// cells_src and cells_tgt will (even after cloning)
  /// refer to different cell storages!
  static void clone(const MeshBoundarySet &boundary_src, MeshBoundarySet &boundary_tgt);

  /// Reset all data of given dimension
  void clear(const Uint dim);

  /// Create one domain domains
  bdry_facets_shared_ptr create(const Uint dim, const std::string &domain_name);

  /// Get one domain by name
  bdry_facets_shared_ptr domain(const std::string &domain_name);

  /// Get one domain by name, const version
  const bdry_facets_shared_ptr domain(const std::string &domain_name) const;

  /// Get one domain by id and by dimension
  bdry_facets_shared_ptr domain(const std::string &domain_name, const Uint dim);

  /// Get one domain by id and dimension, const version
  const bdry_facets_shared_ptr domain(const std::string &domain_name, const Uint dim) const;

  /// Return an iterator range covering one part of boundary
  const common::IteratorRange<const_bdry_topo_iterator> domain_range(
      const std::string &domain_name) const;

  /// Return an iterator range covering one part of boundary
  const common::IteratorRange<const_bdry_topo_iterator> domain_range(const std::string &domain_name,
                                                                     const Uint dim) const;

  /// Return the total number of domains
  Uint nb_domains() const;

  /// Return the total number of cells in all domains
  Uint total_nb_cells() const;

  /// Return a vector of pointers to all domains
  const std::vector<bdry_facets_shared_ptr> &all_domains() const;

  /// List all domains, this is for debugging
  void list() const;

  private:
  /// Pointer to cell data storage
  common::PtrHandle<mesh_cells_t const> m_mesh_cells;

  /// Vector storing all domains
  std::vector<bdry_facets_shared_ptr> m_domains;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshBoundarySet<MeshConfig>::MeshBoundarySet()
{
  m_domains.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
MeshBoundarySet<MeshConfig>::~MeshBoundarySet()
{
  clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshBoundarySet<MeshConfig>::init(
    typename MeshBoundarySet<MeshConfig>::mesh_cells_t const &mesh_cells)
{
  m_mesh_cells = common::PtrHandle<mesh_cells_t const>(&mesh_cells);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshBoundarySet<MeshConfig>::clear()
{
  //  for(Uint i = 0; i < m_domains.size(); ++i)
  //  {
  //    m_domains[i]->clear();
  //  }
  m_mesh_cells.reset(nullptr);
  m_domains.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshBoundarySet<MeshConfig>::clone(const MeshBoundarySet &boundary_src,
                                        MeshBoundarySet &boundary_tgt)
{
  /*
  using cells_type = typename result_of::cells<MeshConfig>::type;

  boundary_tgt.m_domains.clear();

  for (Uint i = 0; i < boundary_src.m_domains.size(); ++i)
  {
    cells_shared_ptr domain_src = boundary_src.m_domains[i];
    cells_shared_ptr domain_tgt(new cells_type());

    domain_tgt->init(*boundary_tgt.m_cells, domain_src->dim(),
  domain_src->name()); cells_type::clone(*domain_src, *domain_tgt);

    boundary_tgt.m_domains.push_back(domain_tgt);
  }
  */
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshBoundarySet<MeshConfig>::clear(const Uint dim)
{
  std::vector<bdry_facets_shared_ptr> remaining_domains;

  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if (m_domains[i]->dim() != dim)
    {
      remaining_domains.push_back(m_domains[i]);
    }
  }

  m_domains.swap(remaining_domains);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr MeshBoundarySet<MeshConfig>::create(
    const Uint dim, const std::string &domain_name)
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if ((m_domains[i]->dim() == dim) && (m_domains[i]->name() == domain_name))
    {
      std::cerr << "SubCells::create:: domain [\"" << domain_name << ", dim = " << dim
                << "] already exists" << std::endl;

      bdry_facets_shared_ptr dom_ptr;
      return dom_ptr;
    }
  }

  bdry_facets_shared_ptr new_domain(
      new BoundaryFacets<MeshConfig, MeshConfig::TDIM - 1>(m_mesh_cells, domain_name));

  // new_domain->init(*m_cells, dim, domain_name);
  // new_domain->init(dim, domain_name);

  m_domains.push_back(new_domain);
  return new_domain;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr MeshBoundarySet<MeshConfig>::domain(
    const std::string &domain_name)
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if (m_domains[i]->name() == domain_name)
    {
      return m_domains[i];
    }
  }

  bdry_facets_shared_ptr empty_domain;
  return empty_domain;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr const MeshBoundarySet<
    MeshConfig>::domain(const std::string &domain_name) const
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if (m_domains[i]->name() == domain_name)
    {
      return m_domains[i];
    }
  }

  bdry_facets_shared_ptr empty_domain;
  return empty_domain;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr MeshBoundarySet<MeshConfig>::domain(
    const std::string &domain_name, const Uint dim)
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if ((m_domains[i]->name() == domain_name) && (m_domains[i]->dim() == dim))
    {
      return m_domains[i];
    }
  }

  std::cerr << "Didn't find domain [name = " << domain_name << ", dim = " << dim << "]"
            << std::endl;
  bdry_facets_shared_ptr empty_domain;
  return empty_domain;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr const MeshBoundarySet<
    MeshConfig>::domain(const std::string &domain_name, const Uint dim) const
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if ((m_domains[i]->name() == domain_name) && (m_domains[i]->dim() == dim))
    {
      return m_domains[i];
    }
  }

  std::cerr << "Didn't find domain [id = " << domain_name << ", dim = " << dim << "]" << std::endl;
  bdry_facets_shared_ptr empty_domain;
  return empty_domain;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::IteratorRange<typename MeshBoundarySet<MeshConfig>::const_bdry_topo_iterator>
MeshBoundarySet<MeshConfig>::domain_range(const std::string &domain_name) const
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if ((m_domains[i]->name() == domain_name))
    {
      return common::IteratorRange<const_bdry_topo_iterator>(m_domains[i]->cbegin(),
                                                             m_domains[i]->cend());
    }
  }

  std::cerr << "Didn't find domain [id = " << domain_name << "]" << std::endl;
  common::IteratorRange<const_bdry_topo_iterator> empty_range;
  return empty_range;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::IteratorRange<typename MeshBoundarySet<MeshConfig>::const_bdry_topo_iterator>
MeshBoundarySet<MeshConfig>::domain_range(const std::string &domain_name, const Uint dim) const
{
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    if ((m_domains[i]->name() == domain_name) && (m_domains[i]->dim() == dim))
    {
      return common::IteratorRange<const_bdry_topo_iterator>(m_domains[i]->cbegin(*m_mesh_cells),
                                                             m_domains[i]->cend(*m_mesh_cells));
    }
  }
  std::cerr << "Didn't find domain [id = " << domain_name << ", dim = " << dim << "]" << std::endl;
  common::IteratorRange<const_bdry_topo_iterator> empty_range;
  return empty_range;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint MeshBoundarySet<MeshConfig>::nb_domains() const
{
  return m_domains.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint MeshBoundarySet<MeshConfig>::total_nb_cells() const
{
  Uint nb_cells = 0;

  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    nb_cells += m_domains[i]->nb_active_cells();
  }

  return nb_cells;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::vector<typename MeshBoundarySet<MeshConfig>::bdry_facets_shared_ptr> const &MeshBoundarySet<
    MeshConfig>::all_domains() const
{
  return m_domains;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void MeshBoundarySet<MeshConfig>::list() const
{
  std::cout << "Domains (total " << m_domains.size() << " domains): " << std::endl;
  for (Uint i = 0; i < m_domains.size(); ++i)
  {
    std::cout << "  " << m_domains[i]->name() << ", " << m_domains[i]->dim() << "D" << std::endl;
    // std::cout << "Element types in domain " << m_domains[i]->name() <<
    // std::endl; m_domains[i]->print_cell_types(); std::cout << std::endl
    // << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
