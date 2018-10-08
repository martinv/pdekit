#ifndef PDEKIT_Mesh_Containers_Dof_Map_hpp
#define PDEKIT_Mesh_Containers_Dof_Map_hpp

#include <algorithm>
#include <cmath>
#include <ctime>
#include <forward_list>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "common/Component.hpp"
#include "common/IteratorRange.hpp"
#include "mesh/CellBuffer.hpp"
#include "mesh/CellMarker.hpp"
#include "mesh/DofCoordinates.hpp"
#include "mesh/TopologyPredicates.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/adaptation/MeshAdaptSequence.hpp"
#include "mesh/iterators/CellDofIterator.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"
#include "mesh/view/CellDofView.hpp"

namespace pdekit
{

namespace mesh
{

template <typename MeshConfig>
class DofMap : public common::Component
{
  public:
  /// TYPEDEFS:
  using config_t = MeshConfig;

  using shared_ptr = std::shared_ptr<DofMap<MeshConfig>>;

  using dof_iterator       = CellDofIterator<CellDofView<DofMap, ViewIsNotConst>>;
  using const_dof_iterator = CellDofIterator<CellDofView<DofMap, ViewIsConst>>;

  /*
  using dof_iterator_typed = DofIteratorTyped<DofViewTyped<DofMap,
  ViewIsNotConst>>; using const_dof_iterator_typed =
  DofIteratorTyped<DofViewTyped<DofMap, ViewIsConst>>;
  */

  using dof_iterator_typed =
      CellDofIterator<CellDofView<DofMap, ViewIsNotConst>, DofIterFilterTyped>;
  using const_dof_iterator_typed =
      CellDofIterator<CellDofView<DofMap, ViewIsConst>, DofIterFilterTyped>;

  using const_dof_range_typed = common::IteratorRange<const_dof_iterator_typed>;
  using all_dof_type_ranges   = std::vector<const_dof_range_typed>;

  enum
  {
    GDIM = MeshConfig::GDIM
  };

  enum
  {
    TDIM = MeshConfig::TDIM
  };

  public:
  /// Constructor
  DofMap(Tria<MeshConfig> const *tria, const std::string &name);

  /// Initialize the object - set the dimension of entities it holds
  /// and its name
  void init(const Uint dim);

  /// Set the dimension of these cells
  // void init(const Uint dim, const std::string& name);

  /// Copy constructor
  /// @param other cell connectivity from which this cell connectivity
  ///        should be copied
  ///        COPY CONSTRUCTOR IS INTENTIONALLY DISABLED
  DofMap(const DofMap &other) = delete;

  /// Assignment operator (*this) = other
  /// @param  other connectivity from which we are assigning
  /// @return reference to 'this' connectivity
  ///        ASSIGNEMENT OPERATOR IS INTENTIONALLY DISABLED
  DofMap &operator=(const DofMap &other) = delete;

  /// Destructor
  ~DofMap() override = default;

  /// Clear all internal data
  void clear();

  /// Return the dimension of this dof storage
  Uint dim() const;

  /// Return the type name of this class (to use with component)
  std::string derived_type_name() const override;

  /// Return true if the DofMap has been adapted
  bool is_adapted() const;

  /// Return the underlying topology cell
  const CellTopologyView<MeshConfig> tcell(const ActiveIdx idx) const;

  /// Return the underlying topology cell
  const CellTopologyView<MeshConfig> tcell(const FlatIdx idx) const;

  /// Return one active cell on position i
  /// @note index 'i' refers only to position among ACTIVE cells
  const MeshEntity active_cell(const ActiveIdx idx) const;

  /// Return one cell on position i, active or not
  const MeshEntity cell(const FlatIdx idx) const;

  /// Return the type id of i-th cell, selecting only among active cells
  PointSetTag active_cell_std_region_id(const ActiveIdx idx) const;

  /// Fill entity so that it corresponds to the cell on position i
  void fill_cell(const ActiveIdx idx, MeshEntity &entity) const;

  /// Fill entity so that it corresponds to the cell on position i
  void fill_cell(const FlatIdx idx, MeshEntity &entity) const;

  /// Total number of cells in this connectivity
  Uint nb_active_cells() const;

  /// Total number of cells in this connectivity
  Uint nb_all_cells() const;

  /// Total number of nodes in this connectivity
  Uint nb_nodes() const;

  /// Return the number of types that this connectivity holds
  Uint nb_all_cell_types() const;

  /// Clear the existing stored information and re-create new data from
  /// supplied cell buffer. This destroys any previous allocated data!
  void create_from_cells(const CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> &cell_buffer);

  /// Renumber all dofs held in the handler
  template <typename T>
  void renumber_dofs(const std::vector<T> &renumbering);

  /// Renumber all dofs held in the handler by treating cells as blocks
  /// of consecutive dofs
  template <typename T>
  void renumber_dofs_blockwise(const std::vector<T> &cell_renumbering);

  /// Renumber all dofs held in the handler by treating cells as blocks
  /// of consecutive dofs
  /// At the same time, export node reordering that this implies
  template <typename T1, typename T2>
  void renumber_dofs_blockwise(const std::vector<T1> &cell_renumbering,
                               std::vector<T2> &node_renumbering);

  /// Insert new cells into the container
  void adapt(const MeshAdaptSequence<MeshConfig> &adapter);

  /// Update the dof handler following update in the mesh
  void adapt_update(const Tria<MeshConfig> &adapted_mesh);

  /// Update the dof handler following update in the mesh
  /// @param adapted_mesh        ... mesh topology after adaptation
  /// @param active_p_orders     ... polynomial orders of active cells
  void adapt_update(const Tria<MeshConfig> &adapted_mesh,
                    const common::ArrayView<const Uint, _1D, Uint> &active_p_orders);

  /// Return the number of tags that the cells contain
  Uint nb_cell_tag_types() const;

  /// Return the material tag of one cell
  Uint active_cell_tag(const ActiveIdx cell_idx) const;

  /// Return the status of one cell (not necessarily active, going only
  /// over active cells would be useless)
  EntityStatus cell_status(const FlatIdx cell_idx) const;

  /// Get the tag names and values
  void all_tag_names(std::map<Uint, std::string> &names) const;

  /// Set the physical tag of all cells
  void tag_all_active_cells(const std::vector<std::pair<Uint, std::string>> &cell_tag_names,
                            std::vector<Uint> &tag_values);

  /// Get the cell reordering
  const std::vector<Uint> &cell_reordering() const;

  /// Fill an std::vector with all (const) cell groups that
  /// this cell connectivity has
  const all_dof_type_ranges &all_active_dof_groups() const;

  /// Print all cell types stored here:
  void print_cell_types() const;

  /// Create an iterator pointing to the beginning of the cell array
  // type_iterator begin(const Uint cell_type);

  /// Create an iterator pointing to the beginning of the cell array, const
  /// version
  const_dof_iterator_typed begin_typed(const PointSetTag standard_region_tag) const;

  /// Return const iterator regardless of whether Cells is a constant
  /// object or not
  const_dof_iterator_typed cbegin_typed(const PointSetTag standard_region_tag) const;

  /// Create an iterator pointing 1 position after the end of the cell array
  // type_iterator end(const Uint cell_type);

  /// Create an iterator pointing 1 position after the end of the cell array
  const_dof_iterator_typed end_typed(const PointSetTag standard_region_tag) const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_dof_iterator_typed cend_typed(const PointSetTag standard_region_tag) const;

  /// Create an iterator pointing to the beginning of the cell array, const
  /// version
  const_dof_iterator begin() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_dof_iterator cbegin() const;

  /// Create an iterator pointing 1 position after the end of the cell array
  const_dof_iterator end() const;

  /// Return const iterator regardless of whether Cells is a constant object
  /// or not
  const_dof_iterator cend() const;

  common::IteratorRange<const_dof_iterator> as_range() const;

  /// @todo: MAKE THIS METHOD PRIVATE AND LET ONLY FRIENDS USE IT
  /// @todo: DESCRIBE HERE WHAT THE METHOD DOES
  void build_skeleton(const Uint dim_skeleton, std::vector<IncidenceEntry> &facet_id_pairs,
                      std::vector<EntityDofRealign> &facet_permutations,
                      std::vector<Uint> &facet_block_offsets) const;

  /// Upgrade topology to a new polynomial order and reference topology type
  /// in one element
  void upgrade(Tria<MeshConfig> &cell_connectivity, const Uint order,
               const PointSetID ref_topology_type);

  /// Clone a topology into a new one, optionally changing the polynomial
  /// order and local topology type of each cell in target mesh. Both input
  /// and output dof handlers are continuous
  static void clone_continuous(Tria<MeshConfig> &cell_connectivity, DofMap &cell_dofs_in,
                               DofMap &cell_dofs_out, const Uint order,
                               const PointSetID ref_topology_type);

  /// Clone a topology into a new one, optionally changing the polynomial
  /// order and local topology type of each cell in target mesh. The input dof
  /// handler is continuous, the output dof handler is discontinuous
  static void clone_discontinuous(const Tria<MeshConfig> &cell_connectivity, DofMap &cell_dofs_in,
                                  DofMap &cell_dofs_out, const Uint order,
                                  const PointSetID ref_topology_type);

  /// Copy one dof handler into another
  static void make_identical_copy(const DofMap &lhs, DofMap &rhs);

  private:
  /// METHODS

  /// Update iterators after the connectivity changes
  void update_all_active_cell_iterators();

  /// Change the types of standard regions for all cells
  void change_std_region_type(const PointSetID new_std_region_type);

  /// Count all degrees of freedom, each DOF counted only once
  Uint count_unique_dofs();

  /// Copy storage pattern from another DofMap
  /// At the same type, insert new dof ids and element types
  void copy_and_upgrade(const DofMap &other_store, std::vector<Uint> &new_dof_ids,
                        const std::vector<PointSetTag> &new_cell_type_ids);

  static void clone_continuous_data_impl(const Tria<MeshConfig> &cell_connectivity,
                                         const DofMap &cell_dofs_in, const Uint order,
                                         const PointSetID ref_topology_type,
                                         std::vector<Uint> &new_connectivity,
                                         std::vector<PointSetTag> &new_cell_type_ids);

  void p_adapt(const std::vector<Uint> &cell_p_orders);

  void h_adapt(const std::vector<CellTransform> &cell_ops);

  /// DATA
  /// Pointer to underlying triangulation (mesh)
  Tria<MeshConfig> const *m_tria;

  /// Topological dimension of the cells in this connectivity array
  Uint m_dim;

  /// Vector of connectivity data
  /// @note Holds only values referring to ACTIVE cells
  std::vector<Uint> m_node_id;

  /// Vector of offsets delimiting blocks of values
  /// that belong to different cells
  /// @note Holds values referring to ALL cells (including deactivated cells)
  std::vector<Uint> m_offset;

  /// Vector of indexes pointing into the offset array
  /// indicating only active cells
  /// @note Holds only values referring to ACTIVE cells
  std::vector<Uint> m_active_cell_pos;

  /// Vector of element types for each cell
  /// @note Holds values referring to ALL cells (including deactivated cells)
  std::vector<StdRegion> m_cell_type;

  /// Vector of flags to indicate if cell is active
  /// @note Holds values referring to ALL cells (including deactivated cells)
  std::vector<EntityStatus> m_status_flag;

  /// Vector determining cell reordering
  /// @note Holds only values referring to ACTIVE cells
  std::vector<Uint> m_cell_reordering;

  /// Number of nodes in DofMap
  Uint m_nb_unique_dofs;

  /// Ranges for all cell types present in this data container
  /// @note Holds only values referring to ACTIVE cells
  all_dof_type_ranges m_group_iterators;

  /// Physical tag for each cell
  CellMarker m_cell_tag;
};

// ----------------------------------------------------------------------------
// Free functions
// ----------------------------------------------------------------------------

/// Get all standard region tags present in dof handler
template <typename MeshConfig>
std::vector<PointSetTag> all_std_reg_tags(const DofMap<MeshConfig> &dof_handler)
{
  std::vector<PointSetTag> collected_tags;

  for (Uint ac = 0; ac < dof_handler.nb_active_cells(); ++ac)
  {
    const MeshEntity active_cell = dof_handler.active_cell(ActiveIdx(ac));
    const PointSetTag candidate  = active_cell.pt_set_id();
    bool found                   = false;

    for (auto tag : collected_tags)
    {
      if (candidate == tag)
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      collected_tags.push_back(candidate);
    }
  }
  return collected_tags;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::vector<std::tuple<PointSetTag, PointSetTag>> all_std_reg_tags(
    const DofMap<MeshConfig> &dof_handler_a, const DofMap<MeshConfig> &dof_handler_b)
{
  std::vector<std::tuple<PointSetTag, PointSetTag>> collected_tags;

  for (Uint ac = 0; ac < dof_handler_a.nb_active_cells(); ++ac)
  {
    const MeshEntity active_cell_a = dof_handler_a.active_cell(ActiveIdx(ac));
    const MeshEntity active_cell_b = dof_handler_b.active_cell(ActiveIdx(ac));
    const std::tuple<PointSetTag, PointSetTag> candidate(active_cell_a.pt_set_id(),
                                                         active_cell_b.pt_set_id());
    bool found = false;

    for (auto tag : collected_tags)
    {
      if (candidate == tag)
      {
        found = true;
        break;
      }
    }
    if (!found)
    {
      collected_tags.push_back(candidate);
    }
  }
  return collected_tags;
}

// ----------------------------------------------------------------------------
// Methods of DofMap
// ----------------------------------------------------------------------------

template <typename MeshConfig>
DofMap<MeshConfig>::DofMap(Tria<MeshConfig> const *tria, const std::string &name)
    : Component(name), m_tria(tria)
{
  init(_0D);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::init(const Uint dim)
{
  m_dim = dim;
  m_node_id.resize(0);
  m_offset.resize(0);
  m_active_cell_pos.resize(0);
  m_cell_type.resize(0);
  m_status_flag.resize(0);
  m_cell_reordering.resize(0);
  m_nb_unique_dofs = 0;
  m_group_iterators.resize(0);
  m_cell_tag.set_type_name("cell_tag");
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::clear()
{
  m_offset.clear();
  m_active_cell_pos.clear();
  m_cell_type.clear();
  m_status_flag.clear();
  m_cell_reordering.clear();
  m_nb_unique_dofs = 0;
  m_group_iterators.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::dim() const
{
  return m_dim;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::string DofMap<MeshConfig>::derived_type_name() const
{
  return "DofMap";
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline bool DofMap<MeshConfig>::is_adapted() const
{
  return (m_active_cell_pos.size() != m_cell_type.size());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> DofMap<MeshConfig>::tcell(const ActiveIdx idx) const
{
  return m_tria->active_cell(idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> DofMap<MeshConfig>::tcell(const FlatIdx idx) const
{
  return m_tria->cell(idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const MeshEntity DofMap<MeshConfig>::active_cell(const ActiveIdx idx) const
{
  const Uint pos_idx = m_active_cell_pos[idx.id()];
  const common::ArrayView<const Uint, _1D, Uint> dofs(m_node_id.data() + m_offset[pos_idx],
                                                      m_offset[pos_idx + 1] - m_offset[pos_idx]);
  return MeshEntity(dofs, idx.id(), m_cell_type[pos_idx]);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const MeshEntity DofMap<MeshConfig>::cell(const FlatIdx idx) const
{
  const Uint i = idx.id();
  const common::ArrayView<const Uint, _1D, Uint> dofs(m_node_id.data() + m_offset[i],
                                                      m_offset[i + 1] - m_offset[i]);
  return MeshEntity(dofs, i, m_cell_type[i]);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline PointSetTag DofMap<MeshConfig>::active_cell_std_region_id(const ActiveIdx idx) const
{
  return m_cell_type[m_active_cell_pos[idx.id()]].get().pt_set_id();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline void DofMap<MeshConfig>::fill_cell(const ActiveIdx idx, MeshEntity &entity) const
{
  const Uint pos_idx                                    = m_active_cell_pos[idx.id()];
  const common::ArrayView<const Uint, _1D, Uint> limits = common::ArrayView<const Uint, _1D, Uint>(
      m_node_id.data() + m_offset[pos_idx], m_offset[pos_idx + 1] - m_offset[pos_idx]);
  entity.reinit(limits, idx.id(), m_cell_type[pos_idx]);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline void DofMap<MeshConfig>::fill_cell(const FlatIdx idx, MeshEntity &entity) const
{
  const Uint i                                          = idx.id();
  const common::ArrayView<const Uint, _1D, Uint> limits = common::ArrayView<const Uint, _1D, Uint>(
      m_node_id.data() + m_offset[i], m_offset[i + 1] - m_offset[i]);
  entity.reinit(limits, i, m_cell_type[i]);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::nb_active_cells() const
{
  return m_active_cell_pos.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::nb_all_cells() const
{
  if (m_offset.size() <= 1)
    return 0;
  return m_offset.size() - 1;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::nb_nodes() const
{
  return m_nb_unique_dofs;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::nb_all_cell_types() const
{
  std::unordered_set<PointSetTag, PointSetTagHasher> cell_types;

  for (Uint i = 0; i < m_active_cell_pos.size(); ++i)
  {
    cell_types.insert(m_cell_type[m_active_cell_pos[i]].get().pt_set_id());
  }
  return cell_types.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::create_from_cells(
    const CellBuffer<MeshConfig::GDIM, MeshConfig::TDIM> &cell_buffer)
{
  // 1) Fill vertex numbers for all cells (container 'm_node_id')
  // 2) At the same time, start filling offsets
  // 3) During the loop over cells, fill also cell types

  m_node_id.assign(cell_buffer.nb_node_entries(), 0u);
  m_offset.assign(cell_buffer.nb_active_cells() + 1, 0u);
  m_active_cell_pos.assign(cell_buffer.nb_active_cells(), 0u);

  for (Uint ac = 0; ac < m_active_cell_pos.size(); ++ac)
  {
    m_active_cell_pos[ac] = ac;
  }

  m_cell_type.resize(cell_buffer.nb_active_cells());

  Uint current_node_idx = 0u;
  for (Uint c = 0; c < cell_buffer.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cell_buffer.active_cell(ActiveIdx(c));
    // Right now, m_offset[c+1] contains the number of vertices of cell c
    // Offset filling is completed in the loop below
    m_offset[c + 1] = cell.nb_vert();
    m_cell_type[c]  = cell.pt_set_id();
    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      m_node_id[current_node_idx++] = cell.vertex(n);
    }
  }

  // Fix offsets: offset of given cell is equal
  // to the offset of previous cell + number of vertices in current cell
  for (Uint i = 1; i < m_offset.size(); ++i)
  {
    m_offset[i] += m_offset[i - 1];
  }

  m_status_flag.assign(cell_buffer.nb_active_cells(), EntityStatus::Active);
  m_cell_reordering.resize(cell_buffer.nb_active_cells());
  std::iota(m_cell_reordering.begin(), m_cell_reordering.end(), 0);

  // 4) Update all iterators correctly
  update_all_active_cell_iterators();

  m_nb_unique_dofs = count_unique_dofs();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename T>
void DofMap<MeshConfig>::renumber_dofs(const std::vector<T> &renumbering)
{
  const Uint my_nb_nodes = nb_nodes();
  if (renumbering.size() != my_nb_nodes)
  {
    std::cerr << "DofMap::renumber_dofs: renumbering vector has wrong size. "
              << "Expected " << my_nb_nodes << ", got " << renumbering.size() << std::endl;
    return;
  }

  // Update the node ids
  for (Uint i = 0; i < m_node_id.size(); ++i)
  {
    m_node_id[i] = renumbering[m_node_id[i]];
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename T>
void DofMap<MeshConfig>::renumber_dofs_blockwise(const std::vector<T> &cell_renumbering)
{
  const Uint my_nb_cells = nb_active_cells();

  if (cell_renumbering.size() != my_nb_cells)
  {
    std::cerr << "DofMap::renumber_dofs_blockwise: renumbering vector has "
                 "wrong size. "
              << "Expected " << my_nb_cells << ", got " << cell_renumbering.size() << std::endl;
    return;
  }

  std::vector<T> cell_renumbering_inverse;
  cell_renumbering_inverse.resize(cell_renumbering.size());

  for (Uint i = 0; i < cell_renumbering.size(); ++i)
  {
    cell_renumbering_inverse[cell_renumbering[i]] = i;
  }

  // Update the order in which coordinates are stored

  Uint dof_idx = 0;
  for (Uint ac = 0; ac < cell_renumbering_inverse.size(); ++ac)
  {
    const Uint cell_pos         = m_active_cell_pos[cell_renumbering_inverse[ac]];
    const Uint cell_nodes_begin = m_offset[cell_pos];
    const Uint cell_nodes_end   = m_offset[cell_pos + 1];
    for (Uint n = cell_nodes_begin; n < cell_nodes_end; ++n)
    {
      m_node_id[n] = dof_idx;
      dof_idx++;
    }
  }

  m_cell_reordering.resize(cell_renumbering.size());
  std::copy(cell_renumbering.begin(), cell_renumbering.end(), m_cell_reordering.begin());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename T1, typename T2>
void DofMap<MeshConfig>::renumber_dofs_blockwise(const std::vector<T1> &cell_renumbering,
                                                 std::vector<T2> &node_renumbering)
{
  const Uint my_nb_cells = nb_active_cells();

  if (cell_renumbering.size() != my_nb_cells)
  {
    std::cerr << "DofMap::renumber_dofs_blockwise: renumbering vector has "
                 "wrong size. "
              << "Expected " << my_nb_cells << ", got " << cell_renumbering.size() << std::endl;
    return;
  }

  std::vector<T1> cell_renumbering_inverse;
  cell_renumbering_inverse.resize(cell_renumbering.size());

  for (Uint i = 0; i < cell_renumbering.size(); ++i)
  {
    cell_renumbering_inverse[cell_renumbering[i]] = i;
  }

  // Resize the node reordering vector
  node_renumbering.resize(m_node_id.size());

  // First 'simulate' the reordering and store it in 'node_renumbering'
  Uint dof_idx = 0;
  for (Uint ac = 0; ac < cell_renumbering_inverse.size(); ++ac)
  {
    const Uint cell_pos         = m_active_cell_pos[cell_renumbering_inverse[ac]];
    const Uint cell_nodes_begin = m_offset[cell_pos];
    const Uint cell_nodes_end   = m_offset[cell_pos + 1];

    for (Uint n = cell_nodes_begin; n < cell_nodes_end; ++n)
    {
      node_renumbering[m_node_id[n]] = dof_idx++;
    }
  }

  // Now actually renumber nodes
  dof_idx = 0;

  for (Uint ac = 0; ac < cell_renumbering_inverse.size(); ++ac)
  {
    const Uint cell_pos         = m_active_cell_pos[cell_renumbering_inverse[ac]];
    const Uint cell_nodes_begin = m_offset[cell_pos];
    const Uint cell_nodes_end   = m_offset[cell_pos + 1];
    for (Uint n = cell_nodes_begin; n < cell_nodes_end; ++n)
    {
      m_node_id[n] = dof_idx;
      dof_idx++;
    }
  }

  m_cell_reordering.resize(cell_renumbering.size());
  std::copy(cell_renumbering.begin(), cell_renumbering.end(), m_cell_reordering.begin());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::adapt(const MeshAdaptSequence<MeshConfig> &schedule)
{
  if (schedule.adapt_type() == AdaptationType::p)
  {
    p_adapt(schedule.cell_poly_orders());
  }
  else if (schedule.adapt_type() == AdaptationType::h)
  {
    h_adapt(schedule.adapt_ops());
  }
  else if (schedule.adapt_type() == AdaptationType::hp)
  {
    std::cerr << "DofMap::adapt: hp adaptation is not implemented." << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::adapt_update(const Tria<MeshConfig> &adapted_mesh)
{
  if (m_offset.size() <= 1)
  {
    std::cerr << "DofMap::adapt::no cells in container => nothing to "
                 "adapt. Aborting."
              << std::endl;
    return;
  }

  std::cout << "DofMap::adapt::starting h-adaptation of handler [" << name() << "]" << std::endl;

  // const Uint nb_old_cells = m_cell_type.size();
  const Uint tot_nb_new_cells = adapted_mesh.nb_all_cells_in_all_levels();

  Uint nb_active_entries    = 0;
  Uint max_refinement_level = 0;

  // Loop over all cells and estimate storage
  // make sure that each cell has its polynomial order and reference
  // topology type set
  // Cells created due to refinement inherit this information from parents

  std::vector<StdRegion> updated_cell_type;
  updated_cell_type.resize(tot_nb_new_cells);

  std::vector<bool> cell_processed;
  cell_processed.resize(tot_nb_new_cells);
  cell_processed.assign(tot_nb_new_cells, false);

  std::vector<Uint> updated_material_id;
  updated_material_id.resize(tot_nb_new_cells);

  // First go through cells on highest refinement level. These should preserve
  // their information about polynomial degree and reference topology
  // throughout refinements

  for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
  {
    const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(FlatIdx(cell_id));
    max_refinement_level = std::max(max_refinement_level, tcell.refinement_level());

    if (tcell.refinement_level() == 0)
    {
      updated_cell_type[tcell.linear_pos_idx().id()] = m_cell_type[tcell.linear_pos_idx().id()];
      updated_material_id[tcell.linear_pos_idx().id()] =
          m_cell_tag.value(tcell.linear_pos_idx().id());
      cell_processed[tcell.linear_pos_idx().id()] = true;

      if (tcell.status() == EntityStatus::Active)
      {
        nb_active_entries += updated_cell_type[tcell.linear_pos_idx().id()].get().nb_nodes();
      }
    }
  }

  for (Uint l = 1; l <= max_refinement_level; ++l)
  {
    for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
    {
      const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(FlatIdx(cell_id));
      if (tcell.refinement_level() == l)
      {
        const CellTopologyView<MeshConfig> tcell_parent = tcell.parent();
        const StdRegion parent_type                     = tcell_parent.std_region();
        const adapt::CellAdaptOp parent_adapt_op        = tcell_parent.cell_adapt_op();

        const PointSetTag parent_tag =
            updated_cell_type[tcell_parent.linear_pos_idx().id()].get().pt_set_id();
        const Uint parent_p              = parent_tag.poly_order();
        const PointSetID parent_ref_topo = parent_tag.ref_topology();

        const Uint parent_material_id = updated_material_id[tcell_parent.linear_pos_idx().id()];

        const std::vector<CellTopologyView<MeshConfig>> children = tcell_parent.children();

        for (Uint child_id = 0; child_id < children.size(); ++child_id)
        {
          const CellTopologyView<MeshConfig> child_tcell = children[child_id];
          if (!cell_processed[child_tcell.linear_pos_idx().id()])
          {
            const ElemShape child_shape = parent_adapt_op.get().child_elem_shape(child_id);
            const PointSetTag child_tag(child_shape, parent_p, parent_ref_topo);
            const StdRegion child_std_region(child_tag);

            updated_cell_type[child_tcell.linear_pos_idx().id()]   = child_std_region;
            updated_material_id[child_tcell.linear_pos_idx().id()] = parent_material_id;
            cell_processed[child_tcell.linear_pos_idx().id()]      = true;

            if (child_tcell.status() == EntityStatus::Active)
            {
              nb_active_entries +=
                  updated_cell_type[child_tcell.linear_pos_idx().id()].get().nb_nodes();
            }
          }
        } // Loop over children

      } // If the cell is on refinement level l
    }
  } // Loop over levels

  Uint nb_processed_cells = 0;
  for (Uint i = 0; i < cell_processed.size(); ++i)
  {
    if (cell_processed[i])
    {
      nb_processed_cells++;
    }
  }

  /*
  std::cout << "DofMap h-update: processed " << nb_processed_cells << "/"
            << cell_processed.size() << " cells." << std::endl;
  std::cout << "Number of active entires (nodes in dof handler)  " <<
  nb_active_entries
            << std::endl;
  */

  std::vector<Uint> updated_offset;
  updated_offset.resize(tot_nb_new_cells + 1);
  updated_offset.assign(tot_nb_new_cells + 1, 0);

  // Compute new coordinates in the updated dof handler
  // Start filling the new offset array

  m_active_cell_pos.reserve(adapted_mesh.nb_active_cells());
  m_active_cell_pos.resize(0);
  m_status_flag.resize(tot_nb_new_cells);

  std::vector<Real> tmp_work;

  for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
  {
    const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(FlatIdx(cell_id));
    m_status_flag[cell_id]                   = tcell.status();

    if (tcell.status() == EntityStatus::Active)
    {
      updated_offset[cell_id + 1] = updated_cell_type[cell_id].get().nb_nodes();

      m_active_cell_pos.push_back(tcell.linear_pos_idx().id());
    } // If this cell is active
    else
    {
      updated_offset[cell_id + 1] = 0;
    }
  } // Loop over all cells

  // Generate the new offsets
  for (Uint i = 1; i < updated_offset.size(); ++i)
  {
    updated_offset[i] += updated_offset[i - 1];
  }

  m_node_id.resize(nb_active_entries);
  for (Uint n = 0; n < m_node_id.size(); ++n)
  {
    m_node_id[n] = n;
  }

  m_offset.swap(updated_offset);
  // m_active_cell_pos was updated in the loop above
  m_cell_type.swap(updated_cell_type);
  // m_status_flag was updated in the loop above

  m_cell_reordering.resize(m_cell_type.size());
  std::iota(m_cell_reordering.begin(), m_cell_reordering.end(), 0);

  m_nb_unique_dofs = count_unique_dofs();

  std::map<Uint, std::string> marker_names;
  m_cell_tag.all_marker_names(marker_names);

  m_cell_tag.fill(marker_names, updated_material_id);

  update_all_active_cell_iterators();

  std::cout << "DofMap::adapt::finished h-adaptation of handler [" << name() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::adapt_update(
    const Tria<MeshConfig> &adapted_mesh,
    const common::ArrayView<const Uint, _1D, Uint> &active_p_orders)
{
  if (m_offset.size() <= 1)
  {
    std::cerr << "DofMap::adapt::no cells in container => nothing to "
                 "adapt. Aborting."
              << std::endl;
    return;
  }

  std::cout << "DofMap::adapt::starting h-adaptation of handler [" << name() << "]" << std::endl;

  // const Uint nb_old_cells = m_cell_type.size();
  const Uint tot_nb_new_cells = adapted_mesh.nb_all_cells_in_all_levels();

  Uint nb_active_entries    = 0;
  Uint max_refinement_level = 0;

  // Loop over all cells and estimate storage
  // make sure that each cell has its polynomial order and reference
  // topology type set
  // Cells created due to refinement inherit this information from parents

  std::vector<StdRegion> updated_cell_type;
  updated_cell_type.resize(tot_nb_new_cells);

  std::vector<bool> cell_processed;
  cell_processed.resize(tot_nb_new_cells);
  cell_processed.assign(tot_nb_new_cells, false);

  std::vector<Uint> updated_material_id;
  updated_material_id.resize(tot_nb_new_cells);

  // First go through cells on highest refinement level. These should preserve
  // their information about polynomial degree and reference topology
  // throughout refinements

  for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
  {
    const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(cell_id);
    max_refinement_level = std::max(max_refinement_level, tcell.refinement_level());

    if (tcell.refinement_level() == 0)
    {
      if (tcell.status() == EntityStatus::Active)
      {
        // Change the cell tag so that it takes the polynomial order
        // prescribed in active_p_orders
        const PointSetTag old_cell_tag = m_cell_type[tcell.linear_pos_idx()].get().pt_set_id();
        const PointSetTag new_cell_tag(old_cell_tag.elem_shape(),
                                       active_p_orders[tcell.active_idx()],
                                       old_cell_tag.ref_topology());

        updated_cell_type[tcell.linear_pos_idx()].change_type(new_cell_tag);
        /*
        updated_cell_type[tcell.linear_pos_idx()] =
        m_cell_type[tcell.linear_pos_idx()];
        */

        nb_active_entries += updated_cell_type[tcell.linear_pos_idx()].get().nb_nodes();
      }
      else
      {
        updated_cell_type[tcell.linear_pos_idx()] = m_cell_type[tcell.linear_pos_idx()];
      }

      updated_material_id[tcell.linear_pos_idx()] = m_cell_tag.value(tcell.linear_pos_idx());
      cell_processed[tcell.linear_pos_idx()]      = true;
    }
  }

  for (Uint l = 1; l <= max_refinement_level; ++l)
  {
    for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
    {
      const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(cell_id);
      if (tcell.refinement_level() == l)
      {
        const CellTopologyView<MeshConfig> tcell_parent = tcell.parent();
        const StdRegion parent_type                     = tcell_parent.std_region();
        const adapt::CellAdaptOp parent_adapt_op        = tcell_parent.cell_adapt_op();

        const PointSetTag parent_tag =
            updated_cell_type[tcell_parent.linear_pos_idx()].get().pt_set_id();
        const Uint parent_p              = parent_tag.poly_order();
        const PointSetID parent_ref_topo = parent_tag.ref_topology();

        const Uint parent_material_id = updated_material_id[tcell_parent.linear_pos_idx()];

        const std::vector<CellTopologyView<MeshConfig>> children = tcell_parent.children();

        for (Uint child_id = 0; child_id < children.size(); ++child_id)
        {
          const CellTopologyView<MeshConfig> child_tcell = children[child_id];
          if (!cell_processed[child_tcell.linear_pos_idx()])
          {
            const ElemShape child_shape = parent_adapt_op.get().child_elem_shape(child_id);
            const PointSetTag child_tag(child_shape, parent_p, parent_ref_topo);
            const StdRegion child_std_region(child_tag);

            if (child_tcell.status() == EntityStatus::Active)
            {
              // Change the cell tag so that it takes the
              // polynomial order prescribed in active_p_orders
              const PointSetTag new_cell_tag(child_shape, active_p_orders[child_tcell.active_idx()],
                                             parent_ref_topo);

              updated_cell_type[child_tcell.linear_pos_idx()].change_type(new_cell_tag);

              nb_active_entries += updated_cell_type[child_tcell.linear_pos_idx()].get().nb_nodes();
            }
            else
            {
              updated_cell_type[child_tcell.linear_pos_idx()] = child_std_region;
            }

            updated_material_id[child_tcell.linear_pos_idx()] = parent_material_id;
            cell_processed[child_tcell.linear_pos_idx()]      = true;
          }
        } // Loop over children

      } // If the cell is on refinement level l
    }
  } // Loop over levels

  Uint nb_processed_cells = 0;
  for (Uint i = 0; i < cell_processed.size(); ++i)
  {
    if (cell_processed[i])
    {
      nb_processed_cells++;
    }
  }

  std::vector<Uint> updated_offset;
  updated_offset.resize(tot_nb_new_cells + 1);
  updated_offset.assign(tot_nb_new_cells + 1, 0);

  // Compute new coordinates in the updated dof handler
  // Start filling the new offset array

  m_active_cell_pos.reserve(adapted_mesh.nb_active_cells());
  m_active_cell_pos.resize(0);
  m_status_flag.resize(tot_nb_new_cells);

  for (Uint cell_id = 0; cell_id < tot_nb_new_cells; ++cell_id)
  {
    const CellTopologyView<MeshConfig> tcell = adapted_mesh.cell(cell_id);
    m_status_flag[cell_id]                   = tcell.status();

    if (tcell.status() == EntityStatus::Active)
    {
      updated_offset[cell_id + 1] = updated_cell_type[cell_id].get().nb_nodes();

      m_active_cell_pos.push_back(tcell.linear_pos_idx());

    } // If this cell is active
    else
    {
      updated_offset[cell_id + 1] = 0;
    }
  } // Loop over all cells

  // Generate the new offsets
  for (Uint i = 1; i < updated_offset.size(); ++i)
  {
    updated_offset[i] += updated_offset[i - 1];
  }

  m_node_id.resize(nb_active_entries);
  for (Uint n = 0; n < m_node_id.size(); ++n)
  {
    m_node_id[n] = n;
  }

  m_offset.swap(updated_offset);
  // m_active_cell_pos was updated in the loop above
  m_cell_type.swap(updated_cell_type);
  // m_status_flag was updated in the loop above

  std::iota(m_cell_reordering.begin(), m_cell_reordering.end(), 0);

  m_nb_unique_dofs = count_unique_dofs();

  std::map<Uint, std::string> marker_names;
  m_cell_tag.all_marker_names(marker_names);

  m_cell_tag.fill(marker_names, updated_material_id);

  update_all_active_cell_iterators();

  std::cout << "DofMap::adapt::finished h-adaptation of handler [" << name() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::nb_cell_tag_types() const
{
  return m_cell_tag.nb_types();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline Uint DofMap<MeshConfig>::active_cell_tag(const ActiveIdx cell_idx) const
{
  const Uint pos_idx = m_active_cell_pos[cell_idx.id()];
  return m_cell_tag.value(pos_idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
EntityStatus DofMap<MeshConfig>::cell_status(const FlatIdx cell_idx) const
{
  return m_status_flag[cell_idx.id()];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::all_tag_names(std::map<Uint, std::string> &names) const
{
  m_cell_tag.all_marker_names(names);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::tag_all_active_cells(
    const std::vector<std::pair<Uint, std::string>> &cell_tag_names, std::vector<Uint> &tag_values)
{
  m_cell_tag.fill(cell_tag_names, tag_values);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<Uint> &DofMap<MeshConfig>::cell_reordering() const
{
  return m_cell_reordering;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const typename DofMap<MeshConfig>::all_dof_type_ranges &DofMap<MeshConfig>::all_active_dof_groups()
    const
{
  return m_group_iterators;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::print_cell_types() const
{
  std::cout << "Registered elements in dof store " << m_dim << "D:" << std::endl;
  for (Uint g = 0; g < m_group_iterators.size(); ++g)
  {
    std::cout << "\tElement type [" << g << "] : ";
    std::cout << m_group_iterators[g].begin()->cell_type().as_string() << std::endl;
  }
}

// ----------------------------------------------------------------------------

// CellData::type_iterator CellData::begin(const Uint cell_type)
//{
//  for(Uint g = 0; g < m_group_iterators.size(); ++g)
//  {
//    if ( m_group_iterators[g].begin()->type_id() == cell_type )
//    {
//      return m_group_iterators[g].begin();
//    }
//  }
//}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator_typed DofMap<MeshConfig>::begin_typed(
    const PointSetTag cell_type) const
{
  for (Uint g = 0; g < m_group_iterators.size(); ++g)
  {
    if (m_group_iterators[g].begin()->cell_type() == cell_type)
    {
      return m_group_iterators[g].begin();
    }
  }
  return const_dof_iterator_typed(*this, ActiveIdx(m_cell_type.size()), cell_type);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator_typed DofMap<MeshConfig>::cbegin_typed(
    const PointSetTag cell_type) const
{
  for (Uint g = 0; g < m_group_iterators.size(); ++g)
  {
    if (m_group_iterators[g].begin()->cell_type() == cell_type)
    {
      return m_group_iterators[g].begin();
    }
  }
  return const_dof_iterator_typed(*this, ActiveIdx(m_cell_type.size()), cell_type);
}

// ----------------------------------------------------------------------------

// CellData::type_iterator CellData::end(const Uint cell_type)
//{

//}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator_typed DofMap<MeshConfig>::end_typed(
    const PointSetTag cell_type) const
{
  for (Uint g = 0; g < m_group_iterators.size(); ++g)
  {
    if (m_group_iterators[g].end()->cell_type() == cell_type)
    {
      return m_group_iterators[g].end();
    }
  }
  return const_dof_iterator_typed(*this, ActiveIdx(m_cell_type.size()), cell_type);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator_typed DofMap<MeshConfig>::cend_typed(
    const PointSetTag cell_type) const
{
  for (Uint g = 0; g < m_group_iterators.size(); ++g)
  {
    if (m_group_iterators[g].end()->cell_type() == cell_type)
    {
      return m_group_iterators[g].end();
    }
  }
  return const_dof_iterator_typed(*this, ActiveIdx(m_cell_type.size()), cell_type);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator DofMap<MeshConfig>::begin() const
{
  return const_dof_iterator(*this, ActiveIdx(0));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator DofMap<MeshConfig>::cbegin() const
{
  return const_dof_iterator(*this, ActiveIdx(0));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator DofMap<MeshConfig>::end() const
{
  return const_dof_iterator(*this, ActiveIdx(m_active_cell_pos.size()));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
typename DofMap<MeshConfig>::const_dof_iterator DofMap<MeshConfig>::cend() const
{
  return const_dof_iterator(*this, ActiveIdx(m_active_cell_pos.size()));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
common::IteratorRange<typename DofMap<MeshConfig>::const_dof_iterator> DofMap<
    MeshConfig>::as_range() const
{
  return common::IteratorRange<const_dof_iterator>(cbegin(), cend());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::build_skeleton(const Uint dim_skeleton,
                                        std::vector<IncidenceEntry> &facet_id_pairs,
                                        std::vector<EntityDofRealign> &facet_permutations,
                                        std::vector<Uint> &facet_block_offsets) const
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

  if (is_adapted())
  {
    std::cerr << "DofMap::build_skeleteon: can build skeleton only on "
                 "non-refined dofs. Exiting."
              << std::endl;
    return;
  }

  using TP = TopologyPredicates;

  clock_t start, end;
  Real elapsed;

  start = clock();

  std::cout << "Building " << dim_skeleton << "D skeleton in " << m_dim << "D cells" << std::endl;

  Uint max_node_id = 0;

  //   for (typename CellsType::const_cell_iterator cell_it = cells.cbegin();
  //   cell_it != cells.cend();
  //        ++cell_it)
  for (Uint c = 0; c < nb_active_cells(); ++c)
  {
    const MeshEntity cell = this->active_cell(ActiveIdx(c));
    max_node_id           = std::max(TP::max_entity_vertex(cell), max_node_id);
  }

  const Uint nb_nodes = max_node_id + 1; // We add 1 because nodes are numbered from 0

  // typedef std::list<std::pair<Uint,Uint> > subcell_list_type;
  using IndicenceInfoTriple          = std::tuple<Uint, SUint, EntityDofRealign>;
  using subcell_list_type            = std::forward_list<IndicenceInfoTriple>;
  using const_sub_cell_list_iterator = subcell_list_type::const_iterator;

  using subcell_list_list_type = std::forward_list<std::forward_list<IndicenceInfoTriple>>;
  using const_sub_cell_list_list_iterator = subcell_list_list_type::const_iterator;
  using sub_cell_list_list_iterator       = subcell_list_list_type::iterator;

  using hash_vector_type = std::vector<subcell_list_list_type>;

  hash_vector_type tmp_subcells(nb_nodes);

  // Variable to count how many skeleton edges/facets have been found
  // in the whole mesh
  // This variable does NOT count the multiplicity of multi-edges or
  // multi-faces
  Uint nb_skeleton_items = 0;

  // Cound how many sub-cells (pairs [cell_id, local_id]) have been inserted
  // including their multiplicity
  // Example: if edge is shared by cells 104 and 207 (it is, say, the first
  // edge in cell 104 and fourth edge in cell 207), then it will be stored as
  // {[104,1],[207,4]} and counts as one inserted sub-entity but two inserted
  // pairs (entries)
  Uint nb_inserted_entries = 0;

  bool found_matching_subcell = false;

  // Loop over all cell groups in master connectivity.

  EntityDofRealign p;

  // Loop over all cells
  for (Uint c = 0; c < nb_active_cells(); ++c)
  {
    const MeshEntity cell = this->active_cell(ActiveIdx(c));

    // Loop over all sub-cells of given dimension of one cell
    // sub-cells of dimension 1 = edges
    // sub-cells of dimension 2 = faces
    // std::cout << "cell = " << cell << std::endl;
    // std::cout << "nb of subelem = " << cell.nb_sub_elements(dim_skeleton)
    // << std::endl;
    for (Uint isub = 0; isub < cell.nb_sub_elements(dim_skeleton); ++isub)
    {
      MeshEntity sub_cell = cell.sub_entity(dim_skeleton, isub);

      found_matching_subcell = false;

      const Uint min_vertex_in_sub = TP::min_entity_vertex(sub_cell);

      // Loop over a list of entities which could possibly match the
      // entity 'sub_cell'
      subcell_list_list_type &subcell_list_list = tmp_subcells[min_vertex_in_sub];

      // THE ONLY INTERESTING CASE IS WHEN THE SUBCELL LIST HAS LENGTH 1:
      // IT MEANS THAT IT CONTAINS ONE FACET, WHICH HAS NOT BEEN MATCHED
      // TO ANY THER FACET YET IF THE LENGTH IS 2, TWO FACETS HAVE BEEN
      // ALREADY MATCHED AND NO OTHER FACET SHOULD MATCH THEM IF THERE IS
      // NO SUBCELL LIST CORRESPONDING TO NODE INDEX min_vertex_in_sub, WE
      // NEED TO CREATE IT (SEE CASE  'if (!found_matching_subcell) { ...
      // }'

      // if (std::distance(sub_cell_list_iter->cbegin(),
      // sub_cell_list_iter->cend()) == 1)
      if (!subcell_list_list.empty()) // THIS IS SLOWER, BUT CORRECT: IMAGINE THAT
                                      // 'cells' ARE 2D FACETS IN 3D MESH: THE SET OF
                                      // ALL FACETS IS NON-MANIFOLD, HENCE MORE THAN 2
                                      // 'cells' CAN BE INCIDENT TO THE SAME EDGE, AND
                                      // THE LENGTH OF
                                      // (*sub_cell_list_iter) LIST CAN BE BIGGER THAN
                                      // 2!!!
      {
        for (sub_cell_list_list_iterator sub_cell_list_iter = subcell_list_list.begin();
             sub_cell_list_iter != subcell_list_list.end(); ++sub_cell_list_iter)
        {
          {
            const IndicenceInfoTriple &id_tuple_in_list = *(sub_cell_list_iter->cbegin());
            MeshEntity test_matching_sub_cell =
                this->active_cell(ActiveIdx(std::get<0>(id_tuple_in_list)));
            test_matching_sub_cell.local_transform(dim_skeleton, std::get<1>(id_tuple_in_list));

            if (TP::entities_match(test_matching_sub_cell, sub_cell, p))
            {
              found_matching_subcell = true;

              subcell_list_type &multi_entity = *sub_cell_list_iter;
              // The new entity cannot be pushed in front, because
              // next time we attempt a match,  the match would be
              // attemped against this new entity. Imagine this
              // situation: entity_a inserted into empty list with
              // identity permutation
              //            entity_b inserted into the list with
              //            flip permutation (i.e. b has to be
              //            flipped in order to have the same
              //            orientation as a) entity_c has the
              //            same permutation as b in order to
              //            match entity_a, but because we are
              //            comparing against the front of the
              //            list (currently entity_b), the
              //            detected permutation will be
              //            'identity', because b and c are
              //            oriented the same way instead of
              //            'flip'/'rotation', which is the
              //            permutation of c to obtain a.
              //            Therefore we keep the same entity in
              //            front to always compute permutation
              //            against the same data and insertion of
              //            new entities is done __after__ the
              //            first element

              // multi_entity.push_front(std::make_tuple(cell.idx(),
              // isub, p));
              multi_entity.insert_after(multi_entity.begin(), std::make_tuple(cell.idx(), isub, p));
              nb_inserted_entries++;
              break;
            }
            else if (TP::entities_match_reverse(test_matching_sub_cell, sub_cell, p))
            {
              found_matching_subcell = true;

              subcell_list_type &multi_entity = *sub_cell_list_iter;
              // The new entity cannot be pushed in front, because
              // next time we attempt a match,  the match would be
              // attemped against this new entity.

              // multi_entity.push_front(std::make_tuple(cell.idx(),
              // isub, p));
              multi_entity.insert_after(multi_entity.begin(), std::make_tuple(cell.idx(), isub, p));
              nb_inserted_entries++;
              break;
            }

          } // If the current list of lists is not empty

        } // Loop over the list of (possibly) matching sub-cells
      }

      if (!found_matching_subcell)
      {
        // subcell_list.push_back(IncidencePair(cell.idx(),isub)); //
        // USE WITH LIST

        subcell_list_list.push_front(subcell_list_type());

        subcell_list_type &new_id_pair_list = *(subcell_list_list.begin());
        EntityDofRealign p;
        p.change_type(sub_cell.pt_set_id(),
                      EntityRealignCode::identity(sub_cell.pt_set_id().elem_shape()));
        new_id_pair_list.push_front(std::make_tuple(cell.idx(), isub, p)); // USE PUSH_FRONT
                                                                           // WITH FORWARD_LIST
        nb_skeleton_items++; // Count the number of sub-entities (edges,
                             // cells)
        nb_inserted_entries++;
      }

    } // Loop over all sub-entities of 'cell'

  } // iterator loop over all master cells

  // We store sub-cells including their multiplicity
  facet_id_pairs.resize(nb_inserted_entries);

  facet_permutations.resize(nb_inserted_entries);
  // Resize fact_block_offsets to size (nb_skeleton_items+1) and
  // assign 0 to all of its entries
  facet_block_offsets.assign(nb_skeleton_items + 1, 0u);

  nb_skeleton_items   = 0;
  nb_inserted_entries = 0;

  for (Uint i = 0; i < tmp_subcells.size(); i++)
  {
    subcell_list_list_type &subcell_list_list = tmp_subcells[i];

    for (const_sub_cell_list_list_iterator sub_cell_list_iter = subcell_list_list.begin();
         sub_cell_list_iter != subcell_list_list.end(); ++sub_cell_list_iter)
    {
      Uint subcell_multiplicity = 0;

      // Loop over all the pairs [parent_cell,local_id] that express the
      // same subcell The number of these pairs is the 'subcell
      // multiplicity'
      for (const_sub_cell_list_iterator sub_cell_iter = sub_cell_list_iter->cbegin();
           sub_cell_iter != sub_cell_list_iter->cend(); ++sub_cell_iter)
      {
        facet_id_pairs[nb_inserted_entries].cell_idx = std::get<0>(*sub_cell_iter);
        facet_id_pairs[nb_inserted_entries].local_id = std::get<1>(*sub_cell_iter);
        facet_permutations[nb_inserted_entries]      = std::get<2>(*sub_cell_iter);
        nb_inserted_entries++;
        subcell_multiplicity++;
      } // Loop over list of pairs

      // Fill the block offsets
      facet_block_offsets[nb_skeleton_items + 1] =
          facet_block_offsets[nb_skeleton_items] + subcell_multiplicity;
      nb_skeleton_items++;
    } // Loop over list of lists
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "DofMap: building of mesh skeleton of dim " << dim_skeleton << " took " << elapsed
            << " s" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::upgrade(Tria<MeshConfig> &cell_connectivity, const Uint order,
                                 const PointSetID ref_topology_type)
{
  std::cout << "DofMap::upgrade_ -> starting topology upgrade" << std::endl;

  if (cell_connectivity.nb_levels() != 1)
  {
    std::cerr << "DofMap::upgrade: can only update the topology" << std::endl
              << "before the adaptation starts. Exiting." << std::endl;
    return;
  }

  // If the mesh is three-dimensional, build the edge skeleton as well
  if (MeshConfig::TDIM == 3) // Don't use _3D enum here, the compiler would issue a warning
  {
    std::unique_ptr<std::vector<IncidenceEntry>> edge_id_pairs(new std::vector<IncidenceEntry>());
    std::unique_ptr<std::vector<EntityDofRealign>> facet_permutations(
        new std::vector<EntityDofRealign>());
    std::unique_ptr<std::vector<Uint>> facet_block_offsets(new std::vector<Uint>());

    build_skeleton(_1D, *edge_id_pairs, *facet_permutations, *facet_block_offsets);
    cell_connectivity.rebuild_skeleton(_1D, std::move(edge_id_pairs), std::move(facet_permutations),
                                       std::move(facet_block_offsets));
  }

  clock_t start, end;
  Real elapsed;

  start = clock();

  // Map [cell type , nb. of cells of this type]
  std::map<PointSetTag, Uint> cell_type_to_count;
  // Map [cell type (before upgrade), cell_type (after upgrade)]
  std::map<PointSetTag, PointSetTag> old_and_new_ref_topology;

  // We loop overt the cells in cell_storage as they are laid out, without
  // using any special iterator that would pick only cells of certain type,
  // for example!!!
  const Uint nb_cells = m_cell_type.size();

  // Count how many cells of each type are present
  for (Uint c = 0; c < nb_cells; ++c)
  {
    MeshEntity cell = this->active_cell(ActiveIdx(c));
    cell_type_to_count[cell.pt_set_id().store_value()]++;
  }

  std::cout << "DofMap::upgrade: input mesh statistics:" << std::endl;
  std::cout << "\tTopological dimension = " << cell_connectivity.dim() << std::endl;
  if (cell_connectivity.dim() > _1D)
  {
    std::cout << "\tNumber of edges = " << cell_connectivity.active_skeleton_size(_1D) << std::endl;
  }
  if (cell_connectivity.dim() == _3D)
  {
    std::cout << "\tNumber of faces = " << cell_connectivity.active_skeleton_size(_2D) << std::endl;
  }
  std::cout << "\tNumber of cells = " << cell_connectivity.nb_active_cells() << std::endl
            << std::endl;

  std::cout << "DofMap::upgrade: will upgrade the following cells:" << std::endl;
  for (std::map<PointSetTag, Uint>::const_iterator it = cell_type_to_count.begin();
       it != cell_type_to_count.end(); ++it)
  {
    std::cout << "   Type: " << PointSetTag(it->first).as_string() << "   ...  " << it->second
              << " cells" << std::endl;
  }

  ElemShape old_eshape;
  Uint old_rt_order;
  PointSetID old_rt_type;
  StdRegion new_ref_elem;

  // Number of entries in the new cell connectivity (after upgrade)
  Uint connectivity_arr_size = 0;

  // Vector for new cell connectivity
  std::vector<Uint> new_connectivity;
  // Vector for new cell types
  std::vector<PointSetTag> new_cell_type_ids(nb_cells);
  // Vector of offsets in the new connectivity
  std::vector<Uint> new_offset(nb_cells + 1, 0);
  // Correspondence table [old node id, new node id]
  std::vector<Uint> old_to_new_node_id;

  // ************************************************************************
  // Fill the map of corresponding element types before and after upgrade
  // Compute how many connectivity entries there will be in the upgraded mesh
  // ************************************************************************

  for (std::map<PointSetTag, Uint>::const_iterator ct_iter = cell_type_to_count.begin();
       ct_iter != cell_type_to_count.end(); ++ct_iter)
  {
    // PointSetTag::decompose_uid(ct_iter->first, old_eshape, old_rt_order,
    // old_rt_type);
    PointSetTag old_tag = ct_iter->first;
    PointSetTag::decompose_into_fields(old_tag, old_eshape, old_rt_order, old_rt_type);

    // const Uint new_ref_topology =
    // RefTopologyTag::keys_to_uid(old_eshape,order,old_rt_type);
    // const Uint new_ref_topology =
    // PointSetTag::keys_to_uid(old_eshape,order,ref_topology_type);

    const PointSetTag new_ref_topology(old_eshape, order, ref_topology_type);

    old_and_new_ref_topology.insert(std::make_pair(ct_iter->first, new_ref_topology.store_value()));

    new_ref_elem.change_type(new_ref_topology);

    connectivity_arr_size += ct_iter->second * new_ref_elem.get().nb_nodes();
  }

  // std::cout << "The new connectivity for " << nb_cells << " cells will have
  // "
  //          << connectivity_arr_size << " node entries" << std::endl;

  // Fill the array with new cell types.

  for (Uint c = 0; c < new_cell_type_ids.size(); ++c)
  {
    const PointSetTag old_ref_topology = active_cell_std_region_id(ActiveIdx(c));
    std::map<PointSetTag, PointSetTag>::const_iterator new_ct_iter =
        old_and_new_ref_topology.find(old_ref_topology);

    new_cell_type_ids[c] = new_ct_iter->second;
  }

  new_connectivity.resize(connectivity_arr_size);

  // Copy P1 nodes from the old cells to the new connectivity array. The P1
  // nodes will have to be renumbered, but this will be done later
  // Fill the node numbers which would correspond to any other than P1 nodes
  // by INVALID_NODE_ID

  Uint current_node_entry = 0;
  // The largest number of a p1 node in the mesh
  Uint max_p1_node_id     = 0;
  Uint nb_p1_vert_in_cell = 0;

  for (const_dof_iterator dof_iter = this->begin(); dof_iter != this->end(); ++dof_iter)
  {
    const MeshEntity cell = dof_iter->mesh_entity();
    nb_p1_vert_in_cell    = 0;

    // Count p1 vertices in cell and put them into 'new_connectivity' with
    // OLD node numbering
    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (cell.vert_is_p1(n))
      {
        new_connectivity[current_node_entry++] = cell.vertex(n);

        max_p1_node_id = std::max(max_p1_node_id, cell.vertex(n));
        nb_p1_vert_in_cell++;
      }
    }

    // Find reference element corresponding to this cell and set all
    // higher-order nodes of the cell as nodes with invalid node id
    // in 'new_connectivity'
    std::map<PointSetTag, PointSetTag>::const_iterator cell_type_iter =
        old_and_new_ref_topology.find(cell.pt_set_id());

    // Set this reference element type to the 'new' (upgraded) cell type
    new_ref_elem.change_type(cell_type_iter->second);

    const Uint nb_ho_nodes = new_ref_elem.get().nb_nodes() - nb_p1_vert_in_cell;

    for (Uint n_ho = 0; n_ho < nb_ho_nodes; ++n_ho)
    {
      new_connectivity[current_node_entry++] = INVALID_NODE_ID;
    }

    // Store the number of nodes of cell i in new_offset[i+1]
    new_offset[cell.idx() + 1] = new_ref_elem.get().nb_nodes();

  } // Loop over all cells

  // for (Uint i = 0; i < nb_cells; ++i)
  //{
  //  MeshEntity entity(new_connectivity.data() + new_offset[i],
  //  new_offset[i+1]
  //  - new_offset[i],
  //                    i, new_cell_type_ids[i]);
  //  std::cout << "[" << entity.idx() << " ]" << entity << std::endl;
  //}

  // Now we know what's the largest p1 node number in the mesh and can resize
  // the array used for node renumbering
  old_to_new_node_id.resize(max_p1_node_id + 1);
  old_to_new_node_id.assign(max_p1_node_id + 1, INVALID_NODE_ID);

  // Loop again over all cells and mark those nodes that are p1:
  for (const_dof_iterator cell_iter = this->begin(); cell_iter != this->end(); ++cell_iter)
  {
    const MeshEntity cell = cell_iter->mesh_entity();

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (cell.vert_is_p1(n))
      {
        old_to_new_node_id[cell.vertex(n)] = 0;
      }
    }
  } // Loop over all cells

  // Assign the p1 nodes of the old mesh consecutive numbers:
  Uint nb_p1_nodes_in_new_mesh = 0;

  for (Uint i = 0; i < old_to_new_node_id.size(); ++i)
  {
    if (old_to_new_node_id[i] == 0)
    {
      old_to_new_node_id[i] = nb_p1_nodes_in_new_mesh++;
    }
  }

  Uint nb_ho_nodes_in_new_mesh = 0;
  // *****************************************************************************
  // Rewrite the p1 nodes numbers in the connectivity table and
  // give some default numbers to the remaining nodes in the new connectivity
  // *****************************************************************************
  for (Uint i = 0; i < new_connectivity.size(); ++i)
  {
    if (new_connectivity[i] != INVALID_NODE_ID)
    {
      new_connectivity[i] = old_to_new_node_id[new_connectivity[i]];
    }
    else
    {
      new_connectivity[i] = nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh;
      nb_ho_nodes_in_new_mesh++;
    }
  }

  // *****************************************************************************
  // Correct the offset: it is the number of nodes of i-th cell + the sum of
  // the nodes of
  // all previous cells ( i.e . cells with index smaller than i)
  // *****************************************************************************
  for (Uint i = 1; i < new_offset.size(); ++i)
  {
    new_offset[i] += new_offset[i - 1];
  }

  // *****************************************************************************
  // At this point, the p1 nodes should be correctly numbered.
  // All higher-order nodes on edges, facets (and inside elements) have a
  // number filled in the 'new_connectivity' array, but a number of a facet
  // node prescribed by the 'left' adjacent element and 'right' adjacent
  // element is not necessarily the same Now we will loop over facets and
  // assign unique node id to all higher-order facet nodes
  // *****************************************************************************

  for (Uint e = 0; e < cell_connectivity.active_skeleton_size(_1D); ++e)
  {
    TraceIncidences multi_edge = cell_connectivity.active_skeleton_entry(_1D, ActiveIdx(e));

    if (multi_edge.size() > 1)
    {
      // Construct the entities from raw arrays. Pick the first edge as
      // the reference edge:
      const Uint iL = multi_edge.cell_id(0);
      const common::ArrayView<const Uint, _1D, Uint> dofs_L(
          new_connectivity.data() + new_offset[iL], new_offset[iL + 1] - new_offset[iL]);

      MeshEntity left(dofs_L, iL, new_cell_type_ids[iL]);
      left.local_transform(_1D, multi_edge.local_id(0));

      // Loop over all remaining edges in this MultipleIncidence and unify
      // the node numbers
      for (Uint j = 1; j < multi_edge.size(); ++j)
      {

        const Uint iR = multi_edge.cell_id(j);

        const common::ArrayView<const Uint, _1D, Uint> dofs_R(
            new_connectivity.data() + new_offset[iR], new_offset[iR + 1] - new_offset[iR]);

        MeshEntity right(dofs_R, iR, new_cell_type_ids[iR]);
        right.local_transform(_1D, multi_edge.local_id(j));

        EntityDofRealign p(std::make_pair(
            right.pt_set_id(), EntityRealignCode::identity(right.pt_set_id().elem_shape())));
        bool entities_match = true;

        // Check if the entities match as they are (without any rotation
        // or flip)
        for (Uint i = 0; i < right.nb_vert(); ++i)
        {
          if (left.vertex(i) != right.vertex(p.get().vertex(i)) && left.vert_is_p1(i))
            entities_match = false;
        }

        // If entities don't match, check if we could have a match using
        // a flip

        if (!entities_match)
        {
          p.change_type(right.pt_set_id(),
                        EntityRealignCode::single_flip(right.pt_set_id().elem_shape()));
          entities_match = true;

          for (Uint i = 0; i < left.nb_vert(); ++i)
          {
            if ((left.vertex(i) != right.vertex(p.get().vertex(i))) && left.vert_is_p1(i))
              entities_match = false;
          }
        }

        // ********************************************************************
        // If entities match, loop over their (high-order) nodes and
        // unify them
        // ********************************************************************
        // The high-order vertices of the left and right facets do not
        // coincide (have different numbers). Therefore we have to
        // modify the higher-order vertex numbers in the left OR the
        // right facet so that they match
        if (entities_match)
        {
          for (Uint n = 0; n < left.nb_vert(); ++n)
          {
            if (!left.vert_is_p1(n)) // If this vertex is not p1
            {
              const Uint vertex_left  = left.vertex(n);
              const Uint vertex_right = right.vertex(p.get().vertex(n));

              // Unify the vertex number to be the same in the
              // right and left edge. The final vert. number
              // should be the same as seen from the left
              // ('master') edge
              for (Uint j = new_offset[iR]; j < new_offset[iR + 1];
                   ++j) // Here we loop over raw connectivity
                        // entries of the future right cell
              {
                if (new_connectivity[j] == vertex_right)
                {
                  new_connectivity[j] = vertex_left;
                  break;
                }
              }

            } // if vertex is not p1

          } // loop over all vertices of edge

        } // If entities match

      } // loop over all edges in multiple incidence

    } // If this multicell has length > 1

  } // Loop over all edges in the mesh

  // If the mesh has two-dimensional facets (i.e. the mesh itself is 3D, we
  // need to unify the high-order face nodes as well)
  if (cell_connectivity.dim() > _2D)
  {
    for (Uint f = 0; f < cell_connectivity.active_skeleton_size(_2D); ++f)
    {
      TraceIncidences multi_facet = cell_connectivity.active_skeleton_entry(_2D, ActiveIdx(f));

      if (multi_facet.size() == 2)
      {
        // Construct the entities from raw arrays:
        const Uint iL = multi_facet.cell_id(0);
        const common::ArrayView<const Uint, _1D, Uint> dofs_L(
            new_connectivity.data() + new_offset[iL], new_offset[iL + 1] - new_offset[iL]);

        MeshEntity left(dofs_L, iL, new_cell_type_ids[iL]);
        left.local_transform(cell_connectivity.dim() - 1, multi_facet.local_id(0));

        const Uint iR = multi_facet.cell_id(1);
        const common::ArrayView<const Uint, _1D, Uint> dofs_R(
            new_connectivity.data() + new_offset[iR], new_offset[iR + 1] - new_offset[iR]);
        MeshEntity right(dofs_R, iR, new_cell_type_ids[iR]);
        right.local_transform(cell_connectivity.dim() - 1, multi_facet.local_id(1));

        EntityRealignCode permutation_code;
        permutation_code.add_flip();
        EntityDofRealign p(std::make_pair(right.pt_set_id(), permutation_code));

        bool entities_match = true;

        // Check if the entities match as they are (without any
        // rotation, using just a flip)

        for (Uint i = 0; i < left.nb_vert(); ++i)
        {
          if ((left.vertex(i) != right.vertex(p.get().vertex(i))) && left.vert_is_p1(i))
            entities_match = false;
        }

        // Now we will rotate the right entity trying to match it with
        // the left one The reference vertex is the first permuted
        // vertex after a plain flip(no rotations) We will repeat
        // rotations until the first permuted vertex is equal to this
        // reference vertex again

        if (entities_match == false)
        {
          const Uint ref_vertex_right = right.vertex(p.get().vertex(0));
          permutation_code.add_rotation();
          p.change_type(right.pt_set_id(), permutation_code);

          while (!entities_match && (ref_vertex_right != right.vertex(p.get().vertex(0))))
          {
            entities_match = true;
            for (Uint n = 0; n < left.nb_vert(); ++n)
            {
              if (left.vert_is_p1(n) && (left.vertex(n) != right.vertex(p.get().vertex(n))))
                entities_match = false;
            }

            if (!entities_match)
            {
              permutation_code.add_rotation();
              p.change_type(right.pt_set_id(), permutation_code);
            }
          }
        }

        // ********************************************************************
        // If entities match, loop over their (high-order) nodes and
        // unify them
        // ********************************************************************
        // The high-order vertices of the left and right facets do not
        // coincide (have different numbers). Therefore we have to
        // modify the higher-order vertex numbers in the left OR the
        // right facet so that they match
        if (entities_match)
        {
          for (Uint n = 0; n < left.nb_vert(); ++n)
          {
            if (!left.vert_is_p1(n)) // If this vertex is not p1
            {
              const Uint vertex_left  = left.vertex(n);
              const Uint vertex_right = right.vertex(p.get().vertex(n));

              // If the vertex number seen from LEFT facet is
              // smaller, updated this vertex in the right cell so
              // that the numbers coincide.
              if (vertex_left < vertex_right)
              {
                for (Uint j = new_offset[iR]; j < new_offset[iR + 1];
                     ++j) // Here we loop over raw connectivity
                          // entries of the future right cell
                {
                  if (new_connectivity[j] == vertex_right)
                  {
                    new_connectivity[j] = vertex_left;
                    break;
                  }
                }
              } // if vertex_left < vertex_right

              // If the vertex number seen from the RIGHT is
              // smaller, update the left cell
              else
              {
                for (Uint j = new_offset[iL]; j < new_offset[iL + 1]; ++j)
                {
                  if (new_connectivity[j] == vertex_left)
                  {
                    new_connectivity[j] = vertex_right;
                    break;
                  }
                }
              } // else update the left cell
            }   // if vertex is not p1
          }     // loop over all vertices of facet

        } // If entities match

      } // If this multicell has length 2

    } // Loop over all facets in the mesh

  } // If TopologyType::TDIM > 2

  // ******************************************************************************
  // We need to ensure that after we made the high-order nodes unique, their
  // numbering is continuous starting from nb_p1_nodes_in_new_mesh to
  // nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh
  // The numbering should be as follows:
  // [0  ...   nb_p1_nodes_in_new_mesh-1] [nb_p1_nodes_in_new_mesh ...
  // (nb_p1_nodes_in_new_mesh+nb_ho_nodes_in_new_mesh-1)]
  // [ --- p1 nodes of all elements --- ] [ --- high-order nodes of all
  // elements
  // --- ]
  // ******************************************************************************

  Uint max_node_id_in_new_mesh = 0;
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    max_node_id_in_new_mesh = std::max(max_node_id_in_new_mesh, new_connectivity[n]);
  }

  // Use an array for renumbering of the high-order nodes
  // To save space, new number of high-order node N will be stored in this
  // array on position N-nb_p1_nodes_in_new_mesh, because we know we will not
  // need to renumber the first nb_p1_nodes_in_new_mesh nodes
  std::vector<Uint> old_to_new_ho_node_id;
  old_to_new_ho_node_id.resize(max_node_id_in_new_mesh + 1 - nb_p1_nodes_in_new_mesh);
  old_to_new_ho_node_id.assign(old_to_new_ho_node_id.size(), INVALID_NODE_ID);

  // ***********************************************
  // Generate the new numbering for high-order nodes
  // ***********************************************
  nb_ho_nodes_in_new_mesh = 0;
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    // We can easily detect a high-order node in the new connectivity: its
    // number must be at least equal to nb_p1_nodes_in_new_mesh
    if (new_connectivity[n] >= nb_p1_nodes_in_new_mesh)
    {
      if (old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh] == INVALID_NODE_ID)
      {
        old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh] =
            nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh;
        nb_ho_nodes_in_new_mesh++;
      }
    }
  }

  // ******************************************************************
  // Assign the correct numbers to high-order nodes in new_connectivity
  // ******************************************************************
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    if (new_connectivity[n] >= nb_p1_nodes_in_new_mesh)
    {
      new_connectivity[n] = old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh];
    }
  }

  // *********************************************************************
  // Put the new cell connectivity and cell types into the cells container
  // *********************************************************************
  m_node_id.clear();
  m_node_id.swap(new_connectivity);

  for (Uint i = 0; i < m_cell_type.size(); ++i)
  {
    m_cell_type[i].change_type(new_cell_type_ids[i]);
  }

  // m_active_flag remains untouched!

  // m_offset should already have correct size
  for (Uint i = 0; i < m_offset.size(); ++i)
  {
    m_offset[i] = 0;
  }

  for (Uint i = 1; i < m_offset.size(); ++i)
  {
    m_offset[i] = m_offset[i - 1] + m_cell_type[i - 1].get().nb_nodes();
  }

  // The 'active cell' array does not change - all cells are still active
  // m_active_cell

  // Update the iterators marking beginnings and ends of all cell groups
  update_all_active_cell_iterators();

  // *********************************************************************
  // Upgrade the cell coordinates
  // *********************************************************************
  m_nb_unique_dofs = count_unique_dofs();

  // **************************************************************************

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "DofMap::upgrade to P" << order << " took " << elapsed << " s" << std::endl;
  std::cout << "DofMap::upgrade -> finished upgrade" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::clone_continuous(Tria<MeshConfig> &cell_connectivity, DofMap &cell_dofs_in,
                                          DofMap &cell_dofs_out, const Uint order,
                                          const PointSetID ref_topology_type)
{
  std::cout << "Starting cloning dof handler {continuous} ... " << std::endl;

  if (cell_connectivity.nb_levels() != 1)
  {
    std::cerr << "DofMap::clone_continuous: can only clone" << std::endl
              << "before the adaptation starts. Exiting." << std::endl;
    return;
  }

  // if (MeshConfig::TDIM == 3) // Don't use _3D enum here, the compiler would
  // issue a warning
  if ((MeshConfig::TDIM == 3) && (cell_connectivity.active_skeleton_size(_1D) == 0))
  {
    std::cout << "DofMap::clone_continuous: I need 1D skeleton to perform "
                 "the cloning."
              << std::endl;

    std::unique_ptr<std::vector<IncidenceEntry>> edge_id_pairs(new std::vector<IncidenceEntry>());
    std::unique_ptr<std::vector<EntityDofRealign>> facet_permutations(
        new std::vector<EntityDofRealign>());
    std::unique_ptr<std::vector<Uint>> facet_block_offsets(new std::vector<Uint>());

    cell_dofs_in.build_skeleton(_1D, *edge_id_pairs, *facet_permutations, *facet_block_offsets);
    cell_connectivity.rebuild_skeleton(_1D, std::move(edge_id_pairs), std::move(facet_permutations),
                                       std::move(facet_block_offsets));
  }

  clock_t start, end;
  Real elapsed;

  start = clock();

  // Vector for new cell connectivity
  std::vector<Uint> new_connectivity;
  // Vector for new cell types
  std::vector<PointSetTag> new_cell_type_ids;

  // Compute new topology data
  clone_continuous_data_impl(cell_connectivity, cell_dofs_in, order, ref_topology_type,
                             new_connectivity, new_cell_type_ids);

  // *********************************************************************
  // Put the new cell connectivity and cell types into the cells container
  // *********************************************************************
  cell_dofs_out.copy_and_upgrade(cell_dofs_in, new_connectivity, new_cell_type_ids);

  // *********************************************************************
  // Upgrade the cell coordinates
  // *********************************************************************

  cell_dofs_out.m_nb_unique_dofs = cell_dofs_out.count_unique_dofs();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "DofMap::clone_continuous with order P" << order << " took " << elapsed << " s"
            << std::endl;
  std::cout << "DofMap::clone_continuous -> finished cloning (continuous)" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::clone_discontinuous(const Tria<MeshConfig> &cell_connectivity,
                                             DofMap &cell_dofs_in, DofMap &cell_dofs_out,
                                             const Uint order, const PointSetID ref_topology_type)
{
  std::cout << "Starting cloning dof handler {discontinuous} ... " << std::endl;

  if (cell_connectivity.nb_levels() != 1)
  {
    std::cerr << "DofMap::clone_discontinuous: can only clone the handler" << std::endl
              << "before the adaptation starts. Exiting." << std::endl;
    return;
  }

  clock_t start, end;
  Real elapsed;

  start = clock();

  std::map<PointSetTag, StdRegion> old_to_new_elem;

  Uint nb_discont_nodes = 0;

  for (Uint c = 0; c < cell_dofs_in.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cell_dofs_in.active_cell(ActiveIdx(c));

    std::map<PointSetTag, StdRegion>::const_iterator it = old_to_new_elem.find(cell.pt_set_id());
    if (it == old_to_new_elem.end())
    {
      const PointSetTag old_std_region_tag = cell.pt_set_id();
      const PointSetTag new_std_region_tag(old_std_region_tag.elem_shape(), order,
                                           ref_topology_type);

      StdRegion new_ref_elem;
      new_ref_elem.change_type(new_std_region_tag);

      old_to_new_elem.insert(std::pair<PointSetTag, StdRegion>(old_std_region_tag, new_ref_elem));

      it = old_to_new_elem.find(cell.pt_set_id());
    }

    nb_discont_nodes += it->second.get().nb_nodes();
  }

  std::cout << "Total number of nodes in discontinuous cell connectivity will be "
            << nb_discont_nodes << std::endl;

  std::vector<Uint> new_connectivity(nb_discont_nodes, 0u);
  std::vector<PointSetTag> new_cell_type_ids(cell_dofs_in.nb_active_cells());

  nb_discont_nodes        = 0;
  Uint nb_processed_cells = 0;

  for (Uint c = 0; c < cell_dofs_in.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cell_dofs_in.active_cell(ActiveIdx(c));

    const std::map<PointSetTag, StdRegion>::const_iterator it =
        old_to_new_elem.find(cell.pt_set_id());

    const StdRegion new_std_region = it->second;

    for (Uint v = 0; v < new_std_region.get().nb_nodes(); ++v)
    {
      new_connectivity[nb_discont_nodes] = nb_discont_nodes;
      nb_discont_nodes++;
    }
    new_cell_type_ids[nb_processed_cells] = new_std_region.get().pt_set_id();
    nb_processed_cells++;
  }

  // *********************************************************************
  // Put the new cell connectivity and cell types into the cells container
  // *********************************************************************
  cell_dofs_out.copy_and_upgrade(cell_dofs_in, new_connectivity, new_cell_type_ids);

  // *********************************************************************
  // Upgrade the cell coordinates
  // *********************************************************************
  cell_dofs_out.m_nb_unique_dofs = cell_dofs_out.count_unique_dofs();

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "DofMap::clone_discontinuous with order P" << order << " took " << elapsed << " s"
            << std::endl;
  std::cout << "DofMap::clone_discontinuous -> finished cloning (discontinuous)" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::make_identical_copy(const DofMap &lhs, DofMap &rhs)
{
  rhs.m_dim = lhs.m_dim;

  rhs.m_node_id.resize(lhs.m_node_id.size());
  rhs.m_node_id = lhs.m_node_id;

  rhs.m_offset.resize(lhs.m_offset.size());
  rhs.m_offset = lhs.m_offset;

  rhs.m_active_cell_pos.resize(lhs.m_active_cell_pos.size());
  rhs.m_active_cell_pos = lhs.m_active_cell_pos;

  rhs.m_cell_type.resize(lhs.m_cell_type.size());
  rhs.m_cell_type = lhs.m_cell_type;

  rhs.m_status_flag.resize(lhs.m_status_flag.size());
  rhs.m_status_flag = lhs.m_status_flag;

  rhs.m_cell_reordering.resize(lhs.m_cell_reordering.size());
  rhs.m_cell_reordering = lhs.m_cell_reordering;

  rhs.m_nb_unique_dofs = lhs.m_nb_unique_dofs;

  rhs.update_all_active_cell_iterators();

  CellMarker::copy(lhs.m_cell_tag, rhs.m_cell_tag);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::update_all_active_cell_iterators()
{
  std::unordered_set<PointSetTag, PointSetTagHasher> all_cell_types;

  for (Uint ac = 0; ac < (m_active_cell_pos.size()); ++ac)
  {
    all_cell_types.insert(m_cell_type[m_active_cell_pos[ac]].get().pt_set_id());
  }

  m_group_iterators.clear();

  // Loop over all cell types present in this connectivity and for each type,
  // set an iterator pointing to the
  // first cell of given type and last cell of given type
  for (std::unordered_set<PointSetTag, PointSetTagHasher>::const_iterator type_iter =
           all_cell_types.begin();
       type_iter != all_cell_types.end(); ++type_iter)
  {
    Uint pos_begin = m_active_cell_pos.size();

    // Find the first cell with given type
    for (Uint ac = 0; ac < m_active_cell_pos.size(); ++ac)
    {
      if (m_cell_type[m_active_cell_pos[ac]].get().pt_set_id() == *type_iter)
      {
        pos_begin = ac;
        break;
      }
    }

    // The last cell with given type is simply one position after the end of
    // the cell data array (like in STL)
    const Uint pos_end = m_active_cell_pos.size();

    // Insert a pair of iterators pointing to the first and last positions
    // Note that the positions with which we construct the dof iterators
    // are the positions IN THE ARRAY 'm_active_cell_pos'

    m_group_iterators.push_back(
        const_dof_range_typed(const_dof_iterator_typed(*this, ActiveIdx(pos_begin), *type_iter),
                              const_dof_iterator_typed(*this, ActiveIdx(pos_end), *type_iter)));
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::change_std_region_type(const PointSetID new_std_region_type)
{
  // Uint eshape, std_region_order, std_region_type;

  for (Uint c = 0; c < m_active_cell_pos.size(); ++c)
  {
    const Uint pos                       = m_active_cell_pos[c];
    const PointSetTag old_std_region_tag = m_cell_type[pos].get().type_id();
    const PointSetTag new_std_region_tag(old_std_region_tag.elem_shape(),
                                         old_std_region_tag.poly_order(), new_std_region_type);

    m_cell_type[pos].change_type(new_std_region_tag);
  }

  update_all_active_cell_iterators();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint DofMap<MeshConfig>::count_unique_dofs()
{
  std::unordered_set<Uint> unique_dofs;
  for (const auto dof_id : m_node_id)
  {
    unique_dofs.insert(dof_id);
  }
  return unique_dofs.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::copy_and_upgrade(const DofMap &other_store, std::vector<Uint> &new_dof_ids,
                                          const std::vector<PointSetTag> &new_cell_type_ids)
{
  // Check that the new_cell_type_ids have the same size as the current
  // cell_type_ids:
  if (new_cell_type_ids.size() != other_store.m_cell_type.size())
  {
    std::cerr << "DofMap::copy_and_upgrade: the supplied vector of new "
                 "cell types"
              << " has different length than is the number of cells in the "
                 "dof store "
                 "from which"
                 " the storage pattern should be copied. "
              << "Aborting." << std::endl;
    return;
  }

  m_dim = other_store.m_dim;

  m_node_id.clear();
  m_node_id.swap(new_dof_ids);

  m_cell_type.resize(other_store.m_cell_type.size());
  for (Uint i = 0; i < other_store.m_cell_type.size(); ++i)
  {
    m_cell_type[i].change_type(new_cell_type_ids[i]);
  }

  m_offset.resize(other_store.m_offset.size());
  m_offset[0] = 0;

  for (Uint i = 1; i < m_offset.size(); ++i)
  {
    m_offset[i] = m_offset[i - 1] + m_cell_type[i - 1].get().nb_nodes();
  }

  m_active_cell_pos.resize(other_store.m_active_cell_pos.size());
  std::copy(other_store.m_active_cell_pos.begin(), other_store.m_active_cell_pos.end(),
            m_active_cell_pos.begin());

  m_status_flag.resize(other_store.m_status_flag.size());
  std::copy(other_store.m_status_flag.begin(), other_store.m_status_flag.end(),
            m_status_flag.begin());

  m_cell_reordering.resize(other_store.m_cell_reordering.size());
  std::copy(other_store.m_cell_reordering.begin(), other_store.m_cell_reordering.end(),
            m_cell_reordering.begin());

  CellMarker::copy(other_store.m_cell_tag, m_cell_tag);

  update_all_active_cell_iterators();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::clone_continuous_data_impl(const Tria<MeshConfig> &cell_connectivity,
                                                    const DofMap &cell_dofs_in, const Uint order,
                                                    const PointSetID ref_topology_type,
                                                    std::vector<Uint> &new_connectivity,
                                                    std::vector<PointSetTag> &new_cell_type_ids)
{
  clock_t start, end;
  Real elapsed;

  start = clock();

  // Map [cell type , nb. of cells of this type]
  std::map<PointSetTag, Uint> cell_type_to_count;
  // Map [cell type (before upgrade), cell_type (after upgrade)]
  std::map<PointSetTag, PointSetTag> old_and_new_ref_topology;

  // We loop overt the cells in cell_storage as they are laid out, without
  // using any special iterator that would pick only cells of certain type,
  // for example!!!
  const Uint nb_cells = cell_dofs_in.nb_active_cells();

  // Count how many cells of each type are present
  for (Uint c = 0; c < nb_cells; ++c)
  {
    MeshEntity cell = cell_dofs_in.active_cell(ActiveIdx(c));
    cell_type_to_count[cell.pt_set_id().store_value()]++;
  }

  std::cout << "DofMap::input mesh statistics:" << std::endl;
  std::cout << "\tTopological dimension = " << cell_connectivity.dim() << std::endl;
  if (cell_connectivity.dim() > _1D)
  {
    std::cout << "\tNumber of edges = " << cell_connectivity.active_skeleton_size(_1D) << std::endl;
  }
  if (cell_connectivity.dim() == _3D)
  {
    std::cout << "\tNumber of faces = " << cell_connectivity.active_skeleton_size(_2D) << std::endl;
  }
  std::cout << "\tNumber of cells = " << cell_connectivity.nb_active_cells() << std::endl
            << std::endl;

  std::cout << "DofMap::will upgrade the following cells:" << std::endl;
  for (std::map<PointSetTag, Uint>::const_iterator it = cell_type_to_count.begin();
       it != cell_type_to_count.end(); ++it)
  {
    std::cout << "   Type: " << PointSetTag(it->first).as_string() << "   ...  " << it->second
              << " cells" << std::endl;
  }

  ElemShape old_eshape;
  Uint old_rt_order;
  PointSetID old_rt_type;
  StdRegion new_ref_elem;

  // Number of entries in the new cell connectivity (after upgrade)
  Uint connectivity_arr_size = 0;

  new_cell_type_ids.resize(nb_cells);
  // Vector of offsets in the new connectivity
  std::vector<Uint> new_offset(nb_cells + 1, 0);
  // Correspondence table [old node id, new node id]
  std::vector<Uint> old_to_new_node_id;

  // ************************************************************************
  // Fill the map of corresponding element types before and after upgrade
  // Compute how many connectivity entries there will be in the upgraded mesh
  // ************************************************************************

  for (std::map<PointSetTag, Uint>::const_iterator ct_iter = cell_type_to_count.begin();
       ct_iter != cell_type_to_count.end(); ++ct_iter)
  {
    PointSetTag old_tag = ct_iter->first;
    PointSetTag::decompose_into_fields(old_tag, old_eshape, old_rt_order, old_rt_type);

    const PointSetTag new_ref_topology(old_eshape, order, ref_topology_type);

    old_and_new_ref_topology.insert(std::make_pair(ct_iter->first, new_ref_topology.store_value()));

    new_ref_elem.change_type(new_ref_topology);

    connectivity_arr_size += ct_iter->second * new_ref_elem.get().nb_nodes();
  }

  // std::cout << "The new connectivity for " << nb_cells << " cells will have
  // "
  //          << connectivity_arr_size << " node entries" << std::endl;

  // Fill the array with new cell types.

  for (Uint c = 0; c < new_cell_type_ids.size(); ++c)
  {
    const PointSetTag old_ref_topology = cell_dofs_in.active_cell_std_region_id(ActiveIdx(c));
    std::map<PointSetTag, PointSetTag>::const_iterator new_ct_iter =
        old_and_new_ref_topology.find(old_ref_topology);

    new_cell_type_ids[c] = new_ct_iter->second;
  }

  new_connectivity.resize(connectivity_arr_size);

  // Copy P1 nodes from the old cells to the new connectivity array. The P1
  // nodes will have to be renumbered, but this will be done later
  // Fill the node numbers which would correspond to any other than P1 nodes
  // by INVALID_NODE_ID

  Uint current_node_entry = 0;
  // The largest number of a p1 node in the mesh
  Uint max_p1_node_id     = 0;
  Uint nb_p1_vert_in_cell = 0;

  for (const_dof_iterator cell_iter = cell_dofs_in.begin(); cell_iter != cell_dofs_in.end();
       ++cell_iter)
  {
    const MeshEntity cell = cell_iter->mesh_entity();
    nb_p1_vert_in_cell    = 0;

    // Count p1 vertices in cell and put them into 'new_connectivity' with
    // OLD node numbering
    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (cell.vert_is_p1(n))
      {
        new_connectivity[current_node_entry++] = cell.vertex(n);

        max_p1_node_id = std::max(max_p1_node_id, cell.vertex(n));
        nb_p1_vert_in_cell++;
      }
    }

    // Find reference element corresponding to this cell and set all
    // higher-order nodes of the cell as nodes with invalid node id
    // in 'new_connectivity'
    std::map<PointSetTag, PointSetTag>::const_iterator cell_type_iter =
        old_and_new_ref_topology.find(cell.pt_set_id());

    // Set this reference element type to the 'new' (upgraded) cell type
    new_ref_elem.change_type(cell_type_iter->second);

    const Uint nb_ho_nodes = new_ref_elem.get().nb_nodes() - nb_p1_vert_in_cell;

    for (Uint n_ho = 0; n_ho < nb_ho_nodes; ++n_ho)
    {
      new_connectivity[current_node_entry++] = INVALID_NODE_ID;
    }

    // Store the number of nodes of cell i in new_offset[i+1]
    new_offset[cell.idx() + 1] = new_ref_elem.get().nb_nodes();

  } // Loop over all cells

  // for (Uint i = 0; i < nb_cells; ++i)
  //{
  //  MeshEntity entity(new_connectivity.data() + new_offset[i],
  //  new_offset[i+1]
  //  - new_offset[i],
  //                    i, new_cell_type_ids[i]);
  //  std::cout << "[" << entity.idx() << " ]" << entity << std::endl;
  //}

  // Now we know what's the largest p1 node number in the mesh and can resize
  // the array used for node renumbering
  old_to_new_node_id.resize(max_p1_node_id + 1);
  old_to_new_node_id.assign(max_p1_node_id + 1, INVALID_NODE_ID);

  // Loop again over all cells and mark those nodes that are p1:
  for (const_dof_iterator cell_iter = cell_dofs_in.begin(); cell_iter != cell_dofs_in.end();
       ++cell_iter)
  {
    const MeshEntity cell = cell_iter->mesh_entity();

    for (Uint n = 0; n < cell.nb_vert(); ++n)
    {
      if (cell.vert_is_p1(n))
      {
        old_to_new_node_id[cell.vertex(n)] = 0;
      }
    }
  } // Loop over all cells

  // Assign the p1 nodes of the old mesh consecutive numbers:
  Uint nb_p1_nodes_in_new_mesh = 0;

  for (Uint i = 0; i < old_to_new_node_id.size(); ++i)
  {
    if (old_to_new_node_id[i] == 0)
    {
      old_to_new_node_id[i] = nb_p1_nodes_in_new_mesh++;
    }
  }

  Uint nb_ho_nodes_in_new_mesh = 0;
  // *****************************************************************************
  // Rewrite the p1 nodes numbers in the connectivity table and
  // give some default numbers to the remaining nodes in the new connectivity
  // *****************************************************************************
  for (Uint i = 0; i < new_connectivity.size(); ++i)
  {
    if (new_connectivity[i] != INVALID_NODE_ID)
    {
      new_connectivity[i] = old_to_new_node_id[new_connectivity[i]];
    }
    else
    {
      new_connectivity[i] = nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh;
      nb_ho_nodes_in_new_mesh++;
    }
  }

  // *****************************************************************************
  // Correct the offset: it is the number of nodes of i-th cell + the sum of
  // the nodes of
  // all previous cells ( i.e . cells with index smaller than i)
  // *****************************************************************************
  for (Uint i = 1; i < new_offset.size(); ++i)
  {
    new_offset[i] += new_offset[i - 1];
  }

  // *****************************************************************************
  // At this point, the p1 nodes should be correctly numbered.
  // All higher-order nodes on edges, facets (and inside elements) have a
  // number filled in the 'new_connectivity' array, but a number of a facet
  // node prescribed by the 'left' adjacent element and 'right' adjacent
  // element is not necessarily the same Now we will loop over facets and
  // assign unique node id to all higher-order facet nodes
  // *****************************************************************************

  for (Uint e = 0; e < cell_connectivity.active_skeleton_size(_1D); ++e)
  {
    TraceIncidences multi_edge = cell_connectivity.active_skeleton_entry(_1D, ActiveIdx(e));

    if (multi_edge.size() > 1)
    {
      // Construct the entities from raw arrays. Pick the first edge as
      // the reference edge:
      const Uint iL = multi_edge.cell_id(0);

      const common::ArrayView<const Uint, _1D, Uint> dofs_L(
          new_connectivity.data() + new_offset[iL], new_offset[iL + 1] - new_offset[iL]);
      MeshEntity left(dofs_L, iL, new_cell_type_ids[iL]);

      left.local_transform(_1D, multi_edge.local_id(0));

      // Loop over all remaining edges in this MultipleIncidence and unify
      // the node numbers
      for (Uint j = 1; j < multi_edge.size(); ++j)
      {
        const Uint iR = multi_edge.cell_id(j);

        const common::ArrayView<const Uint, _1D, Uint> dofs_R(
            new_connectivity.data() + new_offset[iR], new_offset[iR + 1] - new_offset[iR]);
        MeshEntity right(dofs_R, iR, new_cell_type_ids[iR]);
        right.local_transform(_1D, multi_edge.local_id(j));

        EntityDofRealign p(std::make_pair(
            right.pt_set_id(), EntityRealignCode::identity(right.pt_set_id().elem_shape())));
        bool entities_match = true;

        // Check if the entities match as they are (without any rotation
        // or flip)
        for (Uint i = 0; i < right.nb_vert(); ++i)
        {
          if (left.vertex(i) != right.vertex(p.get().vertex(i)) && left.vert_is_p1(i))
            entities_match = false;
        }

        // If entities don't match, check if we could have a match using
        // a flip

        if (!entities_match)
        {
          p.change_type(right.pt_set_id(),
                        EntityRealignCode::single_flip(right.pt_set_id().elem_shape()));
          entities_match = true;

          for (Uint i = 0; i < left.nb_vert(); ++i)
          {
            if ((left.vertex(i) != right.vertex(p.get().vertex(i))) && left.vert_is_p1(i))
              entities_match = false;
          }
        }

        // ********************************************************************
        // If entities match, loop over their (high-order) nodes and
        // unify them
        // ********************************************************************
        // The high-order vertices of the left and right facets do not
        // coincide (have different numbers). Therefore we have to
        // modify the higher-order vertex numbers in the left OR the
        // right facet so that they match
        if (entities_match)
        {
          for (Uint n = 0; n < left.nb_vert(); ++n)
          {
            if (!left.vert_is_p1(n)) // If this vertex is not p1
            {
              const Uint vertex_left  = left.vertex(n);
              const Uint vertex_right = right.vertex(p.get().vertex(n));

              // Unify the vertex number to be the same in the
              // right and left edge. The final vert. number
              // should be the same as seen from the left
              // ('master') edge
              for (Uint j = new_offset[iR]; j < new_offset[iR + 1];
                   ++j) // Here we loop over raw connectivity
                        // entries of the future right cell
              {
                if (new_connectivity[j] == vertex_right)
                {
                  new_connectivity[j] = vertex_left;
                  break;
                }
              }

            } // if vertex is not p1

          } // loop over all vertices of edge

        } // If entities match

      } // loop over all edges in multiple incidence

    } // If this multicell has length > 1

  } // Loop over all edges in the mesh

  // If the mesh has two-dimensional facets (i.e. the mesh itself is 3D, we
  // need to unify the high-order face nodes as well)
  if (cell_connectivity.dim() > _2D)
  {
    for (Uint f = 0; f < cell_connectivity.active_skeleton_size(_2D); ++f)
    {
      TraceIncidences multi_facet = cell_connectivity.active_skeleton_entry(_2D, ActiveIdx(f));

      if (multi_facet.size() == 2)
      {
        // Construct the entities from raw arrays:
        const Uint iL = multi_facet.cell_id(0);

        const common::ArrayView<const Uint, _1D, Uint> dofs_L(
            new_connectivity.data() + new_offset[iL], new_offset[iL + 1] - new_offset[iL]);
        MeshEntity left(dofs_L, iL, new_cell_type_ids[iL]);
        left.local_transform(cell_connectivity.dim() - 1, multi_facet.local_id(0));

        const Uint iR = multi_facet.cell_id(1);

        const common::ArrayView<const Uint, _1D, Uint> dofs_R(
            new_connectivity.data() + new_offset[iR], new_offset[iR + 1] - new_offset[iR]);
        MeshEntity right(dofs_R, iR, new_cell_type_ids[iR]);
        right.local_transform(cell_connectivity.dim() - 1, multi_facet.local_id(1));

        EntityRealignCode permutation_code;
        permutation_code.add_flip();
        EntityDofRealign p(std::make_pair(right.pt_set_id(), permutation_code));

        bool entities_match = true;

        // Check if the entities match as they are (without any
        // rotation, using just a flip)

        for (Uint i = 0; i < left.nb_vert(); ++i)
        {
          if ((left.vertex(i) != right.vertex(p.get().vertex(i))) && left.vert_is_p1(i))
            entities_match = false;
        }

        // Now we will rotate the right entity trying to match it with
        // the left one The reference vertex is the first permuted
        // vertex after a plain flip(no rotations) We will repeat
        // rotations until the first permuted vertex is equal to this
        // reference vertex again

        if (entities_match == false)
        {
          const Uint ref_vertex_right = right.vertex(p.get().vertex(0));
          permutation_code.add_rotation();
          p.change_type(right.pt_set_id(), permutation_code);

          while (!entities_match && (ref_vertex_right != right.vertex(p.get().vertex(0))))
          {
            entities_match = true;
            for (Uint n = 0; n < left.nb_vert(); ++n)
            {
              if (left.vert_is_p1(n) && (left.vertex(n) != right.vertex(p.get().vertex(n))))
                entities_match = false;
            }

            if (!entities_match)
            {
              permutation_code.add_rotation();
              p.change_type(right.pt_set_id(), permutation_code);
            }
          }
        }

        // ********************************************************************
        // If entities match, loop over their (high-order) nodes and
        // unify them
        // ********************************************************************
        // The high-order vertices of the left and right facets do not
        // coincide (have different numbers). Therefore we have to
        // modify the higher-order vertex numbers in the left OR the
        // right facet so that they match
        if (entities_match)
        {
          for (Uint n = 0; n < left.nb_vert(); ++n)
          {
            if (!left.vert_is_p1(n)) // If this vertex is not p1
            {
              const Uint vertex_left  = left.vertex(n);
              const Uint vertex_right = right.vertex(p.get().vertex(n));

              // If the vertex number seen from LEFT facet is
              // smaller, updated this vertex in the right cell so
              // that the numbers coincide.
              if (vertex_left < vertex_right)
              {
                for (Uint j = new_offset[iR]; j < new_offset[iR + 1];
                     ++j) // Here we loop over raw connectivity
                          // entries of the future right cell
                {
                  if (new_connectivity[j] == vertex_right)
                  {
                    new_connectivity[j] = vertex_left;
                    break;
                  }
                }
              } // if vertex_left < vertex_right

              // If the vertex number seen from the RIGHT is
              // smaller, update the left cell
              else
              {
                for (Uint j = new_offset[iL]; j < new_offset[iL + 1]; ++j)
                {
                  if (new_connectivity[j] == vertex_left)
                  {
                    new_connectivity[j] = vertex_right;
                    break;
                  }
                }
              } // else update the left cell
            }   // if vertex is not p1
          }     // loop over all vertices of facet

        } // If entities match

      } // If this multicell has length 2

    } // Loop over all facets in the mesh

  } // If TopologyType::TDIM > 2

  // ******************************************************************************
  // We need to ensure that after we made the high-order nodes unique, their
  // numbering is continuous starting from nb_p1_nodes_in_new_mesh to
  // nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh
  // The numbering should be as follows:
  // [0  ...   nb_p1_nodes_in_new_mesh-1] [nb_p1_nodes_in_new_mesh ...
  // (nb_p1_nodes_in_new_mesh+nb_ho_nodes_in_new_mesh-1)]
  // [ --- p1 nodes of all elements --- ] [ --- high-order nodes of all
  // elements
  // --- ]
  // ******************************************************************************

  Uint max_node_id_in_new_mesh = 0;
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    max_node_id_in_new_mesh = std::max(max_node_id_in_new_mesh, new_connectivity[n]);
  }

  // Use an array for renumbering of the high-order nodes
  // To save space, new number of high-order node N will be stored in this
  // array on position N-nb_p1_nodes_in_new_mesh, because we know we will not
  // need to renumber the first nb_p1_nodes_in_new_mesh nodes
  std::vector<Uint> old_to_new_ho_node_id;
  old_to_new_ho_node_id.resize(max_node_id_in_new_mesh + 1 - nb_p1_nodes_in_new_mesh);
  old_to_new_ho_node_id.assign(old_to_new_ho_node_id.size(), INVALID_NODE_ID);

  // ***********************************************
  // Generate the new numbering for high-order nodes
  // ***********************************************
  nb_ho_nodes_in_new_mesh = 0;
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    // We can easily detect a high-order node in the new connectivity: its
    // number must be at least equal to nb_p1_nodes_in_new_mesh
    if (new_connectivity[n] >= nb_p1_nodes_in_new_mesh)
    {
      if (old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh] == INVALID_NODE_ID)
      {
        old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh] =
            nb_p1_nodes_in_new_mesh + nb_ho_nodes_in_new_mesh;
        nb_ho_nodes_in_new_mesh++;
      }
    }
  }

  // ******************************************************************
  // Assign the correct numbers to high-order nodes in new_connectivity
  // ******************************************************************
  for (Uint n = 0; n < new_connectivity.size(); ++n)
  {
    if (new_connectivity[n] >= nb_p1_nodes_in_new_mesh)
    {
      new_connectivity[n] = old_to_new_ho_node_id[new_connectivity[n] - nb_p1_nodes_in_new_mesh];
    }
  }

  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "DofMap::clone_continuous with order P" << order << " took " << elapsed << " s"
            << std::endl;
  std::cout << "DofMap::clone_continuous -> finished cloning (continuous)" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::p_adapt(const std::vector<Uint> &cell_p_orders)
{
  if (m_offset.size() <= 1)
  {
    std::cerr << "DofMap::adapt::no cells in container => nothing to "
                 "adapt. Aborting."
              << std::endl;
    return;
  }

  std::cout << "DofMap::adapt::starting p-adaptation of handler [" << name() << "]" << std::endl;

  Uint tot_nb_dofs = 0;
  std::vector<Uint> new_node_id(0);
  std::vector<Uint> new_offset(m_offset.size(), 0);
  std::vector<StdRegion> new_cell_type(m_cell_type.size());

  // ----------------------------------------------------------------------------
  // PASS 1:
  // Estimate the needed storage
  // Only the number of nodes in the dof handler can potentially change,
  // the number of cells will remain
  // ----------------------------------------------------------------------------

  for (Uint c = 0; c < cell_p_orders.size(); ++c)
  {
    const Uint flat_idx                 = m_active_cell_pos[c];
    const PointSetTag old_cell_type_tag = m_cell_type[flat_idx].get().pt_set_id();

    const Uint p_new = cell_p_orders[c];

    const PointSetTag new_cell_type_tag(old_cell_type_tag.elem_shape(), p_new,
                                        old_cell_type_tag.ref_topology());

    new_cell_type[flat_idx].change_type(new_cell_type_tag);
    tot_nb_dofs += new_cell_type[flat_idx].get().nb_nodes();
  }

  // --------------------------------------------------------------------------
  // PASS 2: resize the internal data and recompute coordinates
  // --------------------------------------------------------------------------

  // std::cout << "Number of cells: " << m_cell_type.size() << std::endl;
  // std::cout << "Number of new dofs: " << tot_nb_dofs << std::endl;

  new_node_id.resize(tot_nb_dofs);
  m_nb_unique_dofs = tot_nb_dofs;

  tot_nb_dofs = 0;

  for (Uint c = 0; c < new_cell_type.size(); ++c)
  {
    if (m_status_flag[c] == EntityStatus::Active)
    {
      new_offset[c + 1] = new_cell_type[c].get().nb_nodes();

      const PointSetTag old_cell_type_tag = m_cell_type[c].get().pt_set_id();
      const PointSetTag new_cell_type_tag = new_cell_type[c].get().pt_set_id();

      const MeshEntity cell = this->cell(FlatIdx(c));

      if (old_cell_type_tag != new_cell_type_tag)
      {
        tot_nb_dofs += new_cell_type[c].get().nb_nodes();
      }
      else
      {
        tot_nb_dofs += m_cell_type[c].get().nb_nodes();
      }

    } // If cell is active
  }   // Loop over cells where new coordinates are computed

  // Fix offsets: offset of given cell is equal
  // to the offset of previous cell + number of vertices in current cell
  for (Uint i = 1; i < new_offset.size(); ++i)
  {
    new_offset[i] += new_offset[i - 1];
  }

  /// Fill the dof storage
  ///@todo Use smarter ordering of nodes, for example Cuthill-McKee
  for (Uint n = 0; n < new_node_id.size(); ++n)
  {
    new_node_id[n] = n;
  }

  m_node_id.swap(new_node_id);
  m_offset.swap(new_offset);

  // m_active_cell_pos does not change: p-adaptation does not affect status of
  // cells m_active_cell_pos

  m_cell_type.swap(new_cell_type);

  // Geometry has already been updated ...

  update_all_active_cell_iterators();

  std::cout << "DofMap::adapt::finished p-adaptation of handler [" << name() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void DofMap<MeshConfig>::h_adapt(const std::vector<CellTransform> &cell_ops)
{
  if (m_offset.size() <= 1)
  {
    std::cerr << "DofMap::adapt::no cells in container => nothing to "
                 "adapt. Aborting."
              << std::endl;
    return;
  }

  std::cout << "DofMap::adapt::starting h-adaptation of handler [" << name() << "]" << std::endl;

  adapt::CellAdaptOp active_cell_adapt_op;
  StdRegion child_elem_type;

  Uint nb_dofs_before_adapt = 0;
  Uint nb_new_dofs          = 0;
  Uint nb_dofs_to_remove    = 0;

  // This includes all cells, including those that are not active!
  const Uint tot_nb_cells_before_adapt = m_cell_type.size();
  Uint nb_new_cells                    = 0;

  // ----------------------------------------------------------------------------
  // PASS 1: estimate the needed storage after adaptation
  //         This loop should not modify any member data yet, because in case
  //         the required adaptation is declared invalid (at the end
  //         of this loop), we don't want to leave the object in invalid state
  // ----------------------------------------------------------------------------

  // This is a loop over ACTIVE cells only
  for (Uint ac = 0; ac < cell_ops.size(); ++ac)
  {
    const Uint pos = m_active_cell_pos[ac];
    nb_dofs_before_adapt += m_cell_type[pos].get().nb_nodes();

    const bool adapt_op_is_refinement = CellTransformTraits::is_refinement(cell_ops[ac]);

    if (adapt_op_is_refinement)
    {
      const PointSetTag parent_type_tag = m_cell_type[pos].get().pt_set_id();
      nb_dofs_to_remove += m_cell_type[pos].get().nb_nodes();

      active_cell_adapt_op.change_type(parent_type_tag.elem_shape(), cell_ops[ac]);

      const Uint nb_children = active_cell_adapt_op.get().nb_child_elems();
      nb_new_cells += nb_children;

      for (Uint c = 0; c < nb_children; ++c)
      {
        // Make sure that the element shape, reference topology and
        // polynomial order of child and parent are the same

        const PointSetTag child_type_tag(active_cell_adapt_op.get().child_elem_shape(c),
                                         parent_type_tag.poly_order(),
                                         parent_type_tag.ref_topology());
        child_elem_type.change_type(child_type_tag);

        nb_new_dofs += child_elem_type.get().nb_nodes();
      }
    } // If this is a refinement operation
  }

  /*
  std::cout << "Number of new cells: " << nb_new_cells << std::endl;
  std::cout << "Number of dofs before adaptation = " << nb_dofs_before_adapt
  << std::endl; std::cout << "Number of new dofs: " << nb_new_dofs <<
  std::endl; std::cout << "Number of dofs to remove: " << nb_dofs_to_remove <<
  std::endl;
  */

  if (nb_dofs_to_remove > nb_dofs_before_adapt)
  {
    std::cerr << "DofMap::h_adaptation failed: \n"
              << "  > number of dofs scheduled to remove: " << nb_dofs_to_remove << "\n"
              << "    which is larger than\n"
              << "  > number of currently present nodes: " << nb_dofs_before_adapt << "\n"
              << "Exiting." << std::endl;
    return;
  }

  std::vector<StdRegion> tmp_child_cell_types(nb_new_cells);
  std::vector<Uint> tmp_child_cell_material_id(nb_new_cells);
  Uint child_idx = 0;

  // ----------------------------------------------------------------------------
  // PASS 2: - store the element child types
  //         - deactivate cells that are marked for refinement
  // ----------------------------------------------------------------------------

  // This is a loop over ACTIVE cells only
  for (Uint ac = 0; ac < cell_ops.size(); ++ac)
  {
    const bool adapt_op_is_refinement = CellTransformTraits::is_refinement(cell_ops[ac]);

    if (adapt_op_is_refinement)
    {
      const Uint pos     = m_active_cell_pos[ac];
      m_status_flag[pos] = EntityStatus::NotActive;

      const PointSetTag parent_type_tag = m_cell_type[pos].get().pt_set_id();
      active_cell_adapt_op.change_type(parent_type_tag.elem_shape(), cell_ops[ac]);

      const Uint nb_children = active_cell_adapt_op.get().nb_child_elems();

      for (Uint c = 0; c < nb_children; ++c)
      {
        // Make sure that the element shape, reference topology and
        // polynomial order of child and parent are the same
        const PointSetTag child_type_tag(active_cell_adapt_op.get().child_elem_shape(c),
                                         parent_type_tag.poly_order(),
                                         parent_type_tag.ref_topology());

        StdRegion child_elem_type(child_type_tag);
        // child_elem_type.change_type(child_type_tag);

        tmp_child_cell_types[child_idx]       = child_elem_type;
        tmp_child_cell_material_id[child_idx] = m_cell_tag.value(c);
        child_idx++;
      }
    } // If this is a refinement operation
  }

  const Uint new_dof_storage_size = nb_dofs_before_adapt + nb_new_dofs - nb_dofs_to_remove;
  m_nb_unique_dofs                = new_dof_storage_size;

  // ----------------------------------------------------------------------------
  // PASS 3: update connectivity
  // ----------------------------------------------------------------------------

  const Uint tot_nb_cells_after_adapt = tot_nb_cells_before_adapt + nb_new_cells;

  m_node_id.resize(new_dof_storage_size);
  for (Uint n = 0; n < new_dof_storage_size; ++n)
  {
    m_node_id[n] = n;
  }

  m_cell_type.resize(tot_nb_cells_after_adapt);
  m_status_flag.resize(tot_nb_cells_after_adapt);

  for (Uint newc = 0; newc < nb_new_cells; ++newc)
  {
    m_cell_type[tot_nb_cells_before_adapt + newc]   = tmp_child_cell_types[newc];
    m_status_flag[tot_nb_cells_before_adapt + newc] = EntityStatus::Active;
  }

  m_offset.resize(tot_nb_cells_after_adapt + 1);
  m_offset[0] = 0;

  Uint nb_active_cells = 0;

  for (Uint c = 0; c < m_cell_type.size(); ++c)
  {
    if (m_status_flag[c] == EntityStatus::Active)
    {
      m_offset[c + 1] = m_cell_type[c].get().nb_nodes();
      nb_active_cells++;
    }
    else
    {
      m_offset[c + 1] = 0;
    }
  }

  // m_active_cell_pos.resize(nb_active_cells, 0);
  m_active_cell_pos.resize(nb_active_cells);
  m_active_cell_pos.assign(nb_active_cells, 0);

  nb_active_cells = 0;
  for (Uint c = 0; c < (m_offset.size() - 1); ++c)
  {
    m_offset[c + 1] += m_offset[c];
    if (m_status_flag[c] == EntityStatus::Active)
    {
      m_active_cell_pos[nb_active_cells] = c;
      nb_active_cells++;
    }
  }

  // Update the cell material ids
  m_cell_tag.push_back_values(tmp_child_cell_material_id);
  update_all_active_cell_iterators();

  /*
  std::cout << "DofMap::adapt::number of all cells: " << m_offset.size() -
  1
  << std::endl; std::cout << "DofMap::adapt::number of active cells: " <<
  m_active_cell_pos.size()
            << std::endl;
  std::cout << "DofMap::adapt::number of material ids: " <<
  m_cell_tag.nb_values() << std::endl;
  */

  std::cout << "DofMap::adapt::finished h-adaptation of handler [" << name() << "]" << std::endl;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::ostream &send_to_output_stream(std::ostream &os, const DofMap<MeshConfig> &cells)
{
  for (Uint c = 0; c < cells.nb_active_cells(); ++c)
  {
    const MeshEntity cell = cells.active_cell(ActiveIdx(c));
    os << "[" << cell.idx() << "] " << cell << std::endl;
  }
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
