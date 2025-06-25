#ifndef PDEKIT_Mesh_Containers_Tria_Cells_hpp
#define PDEKIT_Mesh_Containers_Tria_Cells_hpp

#include <algorithm>
#include <fstream>
#include <iostream>

#include "common/BlockArray.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseMatView.hpp"
#include "mesh/CellGeometry.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/EntityStatus.hpp"
#include "mesh/MeshIndex.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"
#include "mesh/containers/CellPath.hpp"
#include "mesh/local_topology/TraceIncidences.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"

namespace pdekit
{

namespace mesh
{

// Forward declaration so that we can make TopologyCell and MeshTopology friends
// of TriaCells
template <typename MeshConfig>
class CellTopologyView;

namespace internal
{

template <typename MeshConfig>
class MeshTopologyAdaptAlgorithm;

/**
 *  This class stores cells on one level of adaptive mesh. The cells are defined
 * by the indexes of their faces, which are stored in another data structure.
 *  Because the cells are stored in a linear array, offset vector defining
 * position of data related to one cell must be used. In addition, the following
 * information per cell is stored:
 *   - index of parent (which is stored in CellConnectivity one level higher)
 *   - indexes of children (which are stored in CellConnectivity one level
 * lower)
 */

// Forward declaration to enable the output operator
template <typename MeshConfig>
class TriaCells;

template <typename MeshConfig>
std::ostream &operator<<(std::ostream &os, TriaCells<MeshConfig> const &ml);

template <typename MeshConfig>
class TriaCells
{
  public:
  /// TYPEDEFS
  typedef MeshConfig mesh_config;

  /// Default constructor
  TriaCells();

  /// Copy constructor is disabled, so that objects of this class
  /// cannot be accidentally copied
  TriaCells(const TriaCells &other) = delete;

  /// Default destructor
  ~TriaCells();

  /// Assignement operator is disabled
  TriaCells &operator=(const TriaCells &other) = delete;

  /// Get the list of ids of all faces of this cell
  const common::ArrayView<const Uint, _1D, Uint> cell_facet_ids(const FlatIdx cell_idx) const;

  /// Number of all cells in this CellConnectivity
  /// @note: this returns the number of cells regardless of whether they are
  /// active
  ///         or not!
  Uint nb_all_cells() const;

  /// Number of all ACTIVE cells in this connectivity
  Uint nb_active_cells() const;

  /// Return the number of all levels in these cells
  Uint nb_all_levels() const;

  /// Add a new cell
  /// @param cell_type ... StdRegion type of this cell
  /// @return index of the newly created cell in this CellConnectivity
  Uint add_cell(const StdRegion cell_type, const Uint level, const EntityStatus status,
                const Uint parent_pos_idx, const math::DenseDMat<Real> &cell_coords);

  /// Define the neighbours of new cells
  void add_new_cell_neighbours(const Uint facet_numbering_base,
                               const std::vector<Uint> &new_incident_facet_ids,
                               const std::vector<Uint> &new_facet_offsets);

  /// Return one topology cell
  const CellTopologyView<MeshConfig> cell(const FlatIdx idx) const;

  /// Return one ACTIVE topology cell
  const CellTopologyView<MeshConfig> active_cell(const ActiveIdx idx) const;

  /// Update the status of all cells: those that are leafs (have no children)
  /// will be marked as active, all other cells will be not active
  /// Update the list of active cells
  void deactivate_nonleafs_and_update_active_cell_positions();

  /// Set correctly index of a facet which is a part of cell interface
  // void update_cell_facet_idx(const Uint facet_idx, const LocalIndicences&
  // facet);

  /// Get the index of cell's parent, which is stored in CellConnectivity
  /// above this one
  const CellTopologyView<MeshConfig> cell_parent(CellTopologyView<MeshConfig> const &cell) const;

  /// Set the index of cell's parent
  /// @param cell_idx ... index of cell who's parent is being set
  /// @param parent_idx ... index of the parent
  // void set_parent(const Uint cell_idx, const Uint parent_idx);

  /// Get the child indexes of one cell
  const std::vector<CellTopologyView<MeshConfig>> cell_children(const FlatIdx idx) const;

  /// Generate the path to a cell
  /// @param start_tcell     ... topology cell whose path we're trying to
  /// generate
  /// @param zero_level_cell ... root cell of the generated path
  /// @param path_entries    ... child indices of the path, going from
  /// zero_level_cell
  ///                            to path_entries. The first entry of
  ///                            path_entries corresponds to local child id of
  ///                            cell on level 1
  void path(const CellTopologyView<MeshConfig> &start_tcell,
            CellTopologyView<MeshConfig> &zero_level_cell, std::vector<Uint> &path_entries) const;

  /// Get the leaf element of given path
  /// @param  path  ... path to the leaf
  /// @return topology cell which is the last element in the path
  const CellTopologyView<MeshConfig> path_leaf(const CellPath &path) const;

  /// Add children to cell
  /// @param cell_idx ... index of cell to which we are adding children
  /// @param child_indexes ... vector of indexes of the new children
  /// @param parent status ... new status of parent cell
  void add_child_indexes(const Uint cell_idx, std::vector<Uint> const &child_indexes,
                         const EntityStatus &parent_status);

  /// Disconnect the children of given cell from it
  /// Set a new status of the cell and also the status of the children
  /// @param cell_idx ... index of cell to which we are adding children
  /// @param parent status ... new status of parent cell
  void mark_and_disconnect_children(const FlatIdx cell_idx, const EntityStatus &parent_status,
                                    const EntityStatus &child_status);

  /// Remove all leaf cells with given status
  /// If the cell has a given status but also children (i.e. it is not a leaf
  /// cell), then the cell won't be removed.
  /// @param status ... status of cells considered for removal
  void remove_all_leaf_cells_w_status(const EntityStatus &status,
                                      std::vector<FlatIdx> &new_cell_id);

  /// Fill all the data by providing vector of faces,
  /// and a vector of cell types. Note that the length of the vector
  /// of faces must be such that each cell can retrieve its own faces.
  /// For example, when the vector of cell types contains four types 'P1
  /// triangle', then the length of the vector of faces must be 4 x 3 = 12,
  /// since each triangle has 4 faces
  /// @remark: the input vectors will be swapped, with the member variables,
  /// leaving the input
  ///          vectors empty!
  void emplace_cells(std::vector<StdRegion> &cell_types, std::vector<Uint> &cell_faces,
                     common::BlockArray<Real, Uint> const &cell_coordinates);

  /// Rebuild information about neighbours of each cell
  template <typename TraceTopologyIterator>
  void rebuild_neighbour_information(const TraceTopologyIterator &facets_begin,
                                     const TraceTopologyIterator &facets_end);

  /// Change the type of each standard region
  /// @param rule ... a function that takes a PointSetTag as an argument
  ///                 and returns another PointSetTag
  ///                 The input is the existing standard region id
  ///                 The output is the new (desired) standard region tag
  template <typename StdRegChangeRule>
  void change_std_region_types(const StdRegChangeRule &rule);

  /// Perform a geometrical transformation of the mesh
  /// @param rule .. a void function which takes two parameters:
  ///                1) an input vector view (DenseConstVecView) representing
  ///                   the coordinates of one node
  ///                2) an output vector view (DenseConstVecView) representing
  ///                   the new (transformed) coordinates of the same node
  template <typename GeoTransformRule>
  void geo_transform(const GeoTransformRule &rule);

  /// Perform a geometrical transformation of one cell in mesh
  /// @param idx  ... index of cell whose coordinates should be transformed
  /// @param rule .. a void function which takes two parameters:
  ///                1) an input vector view (DenseConstVecView) representing
  ///                   the coordinates of one node
  ///                2) an output vector view (DenseConstVecView) representing
  ///                   the new (transformed) coordinates of the same node
  template <typename GeoTransformRule>
  void geo_transform(const ActiveIdx idx, const GeoTransformRule &rule);

  /// Clear all internal data
  void clear();

  /// Return the position of first active cell
  Uint first_active_cell_pos() const;

  /// Return the position of the last active cell
  Uint last_active_cell_pos() const;

  /// Overloaded output operator to print the contents of the mesh level. This
  /// is used for debugging.
  friend std::ostream &operator<<<MeshConfig>(std::ostream &os, const TriaCells &ml);

  /// Clone one mesh level
  static void clone(const TriaCells &cell_conn_in, TriaCells &cell_conn_out);

  void print_cell_coords_to_file(const std::string &filename) const;

  private:
  /// FRIENDS
  /// This means that TopologyCell is only a friend of TriaCells
  /// when its template parameter matches the template parameter of TriaCells!
  friend class CellTopologyView<MeshConfig>;

  friend class pdekit::mesh::internal::MeshTopologyAdaptAlgorithm<MeshConfig>;

  /// METHODS:

  /// Return the indexes of all faces of one cell, immutable version
  const common::ArrayView<const Uint, _1D, Uint> cell_faces(const FlatIdx cell_idx) const;

  /// Reorder the cell faces
  void reorder_cell_facets(const CellTopologyView<MeshConfig> &cell,
                           const std::vector<Uint> &reordering);

  /// Return the type of given cell
  StdRegion cell_type(const FlatIdx cell_idx) const;

  /// Return adaptation operation that was applied to given cell
  const adapt::CellAdaptOp cell_adapt_op(const FlatIdx cell_idx) const;

  /// Return active index of cell given its flat index
  const ActiveIdx active_idx(const FlatIdx idx) const;

  /// Return the number of children of this cell
  Uint nb_children(const FlatIdx idx) const;

  /// Return cell coordinates of one cell
  const mesh::CellGeometry<MeshConfig::GDIM> cell_geometry(const FlatIdx idx) const;

  /// Return cell coordinates of one cell
  const mesh::CellGeometry<MeshConfig::GDIM> cell_geometry(const FlatIdx idx,
                                                           const Uint sub_entity_dim,
                                                           const Uint local_id) const;

  /// Set the adaptation operation that was applied to given cell
  void set_cell_adapt_op(const FlatIdx cell_idx, const adapt::CellAdaptOp &adapt_op);

  /// Return the level of i-th cell
  /// @param ... index of cell whose level should be returned
  Uint refinement_level(const FlatIdx cell_idx) const;

  /// Return the status of i-th cell
  /// @param ... index of the cell whose status should be returned
  EntityStatus cell_status(const FlatIdx cell_idx) const;

  /// Return the status of all cells
  const std::vector<EntityStatus> &cell_status() const;

  /// Change the status of one cell
  /// @param linear_pos_idx ... linear index of cell whose status we're
  /// changing
  void change_status(const FlatIdx linear_pos_idx, const EntityStatus status);

  /// Check that all cells incident to given 'LocalIncidences' block, i.e. all
  /// cells whose ids are listed in LocalIncidences, have given status
  bool all_incident_cells_have_status(const TraceIncidences &incidences,
                                      const EntityStatus status) const;

  /// Number of all levels contained in this cell connectivity
  Uint m_max_level;

  /// Vector storing the mesh topology by storing the face indexes
  /// of each cell
  std::vector<Uint> m_cell_faces;

  /// Offsets delimiting faces of each cell
  std::vector<Uint> m_cell_face_offset;

  /// Vector of cell types
  std::vector<StdRegion> m_cell_type;

  /// Vector storing adaptive operation that has been applied to each cell
  std::vector<adapt::CellAdaptOp> m_cell_adapt_op;

  /// Vector storing the level for each cell
  std::vector<Uint> m_level;

  /// Vector storing the status of each cell
  std::vector<EntityStatus> m_status;

  /// Position of active cells - this vector holds
  /// indexes of all cells in CellConnectivity (see below)
  /// that are currently marked as active
  /// m_active_to_linear_idx[i] holds the linear (absolute)
  /// position of i-th active cell
  std::vector<Uint> m_active_to_linear_idx;

  /// The inverse mapping of the above: given the
  /// absolute (linear) position, get the index among all active cells
  std::vector<Uint> m_linear_to_active_idx;

  /// Index of parent for each cell. This is index
  /// in the CellConnectivity container one level above.
  std::vector<Uint> m_cell_parent;

  /// List of all children for each cell
  /// @todo: replace this by some more memory-aligned structure
  std::vector<std::vector<Uint>> m_cell_children;

  /// Coordinates of all cells stored in blocks
  /// Each block stores the coordinates of all nodes of each cell
  /// as {(x0,y0,z0), (x1,y1,z1), ..... (xn,yn,zn)}
  common::BlockArray<Real, Uint> m_cell_coords;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
TriaCells<MeshConfig>::TriaCells()
{
  clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
TriaCells<MeshConfig>::~TriaCells()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::ArrayView<const Uint, _1D, Uint> TriaCells<MeshConfig>::cell_facet_ids(
    const FlatIdx cell_idx) const
{
  const Int cell_pos = cell_idx.id();
  const common::ArrayView<const Uint, _1D, Uint> facet_ids(
      m_cell_faces.data() + m_cell_face_offset[cell_pos],
      m_cell_face_offset[cell_pos + 1] - m_cell_face_offset[cell_pos]);
  return facet_ids;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::nb_all_cells() const
{
  return m_cell_type.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::nb_active_cells() const
{
  return m_active_to_linear_idx.size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::nb_all_levels() const
{
  return m_max_level;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::add_cell(const StdRegion cell_type, const Uint level,
                                     const EntityStatus status, const Uint parent_pos_idx,
                                     const math::DenseDMat<Real> &cell_coords)
{
  adapt::CellAdaptOp adapt_op;

  /*
  const Uint nb_faces = cell_type.get().nb_entities(MeshConfig::TDIM - 1);
  for (Uint f = 0; f < nb_faces; ++f)
  {
    m_cell_faces.push_back(INVALID_FACET_ID);
  }
  m_cell_face_offset.push_back(m_cell_face_offset.back() + nb_faces);
  */

  m_cell_type.push_back(cell_type);

  adapt_op.change_type(
      adapt::CellAdaptOpTag(cell_type.get().pt_set_id().elem_shape(), CellTransform::NO_TRANS));
  m_cell_adapt_op.push_back(adapt_op);

  m_level.push_back(level);
  m_status.push_back(status);
  m_cell_parent.push_back(parent_pos_idx);
  m_cell_children.push_back({});

  if (status == EntityStatus::Active)
  {
    m_active_to_linear_idx.push_back(m_cell_type.size() - 1);
    m_linear_to_active_idx.push_back(m_active_to_linear_idx.size() - 1);
  }
  else
  {
    m_linear_to_active_idx.push_back(INVALID_CELL_ID);
  }

  const common::ArrayView<const Real, _1D, Uint> coords_view(
      &cell_coords(0, 0), cell_coords.rows() * cell_coords.cols());

  m_cell_coords.create_back_block(coords_view.size());
  m_cell_coords.fill_last_block(coords_view);

  return m_cell_type.size() - 1;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::add_new_cell_neighbours(const Uint facet_numbering_base,
                                                    const std::vector<Uint> &new_incident_facet_ids,
                                                    const std::vector<Uint> &new_facet_offsets)
{
  if (new_facet_offsets.size() < 2)
  {
    std::cerr << "TriaCells:: add_new_cell_neighbours: no facets to add." << std::endl;
    return;
  }

  const Uint old_cell_faces_size = m_cell_faces.size();
  m_cell_faces.resize(old_cell_faces_size + new_incident_facet_ids.size());

  for (Uint f = 0; f < new_incident_facet_ids.size(); ++f)
  {
    m_cell_faces[old_cell_faces_size + f] = new_incident_facet_ids[f] + facet_numbering_base;
  }

  const Uint old_facet_offsets_size = m_cell_face_offset.size();
  const Uint old_max_offset         = m_cell_face_offset.back();
  m_cell_face_offset.resize(old_facet_offsets_size + new_facet_offsets.size() - 1);

  for (Uint i = 0; i < (new_facet_offsets.size() - 1); ++i)
  {
    m_cell_face_offset[old_facet_offsets_size + i] = old_max_offset + new_facet_offsets[i + 1];
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const CellTopologyView<MeshConfig> TriaCells<MeshConfig>::cell(const FlatIdx idx) const
{
  return CellTopologyView<MeshConfig>(this, idx);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const CellTopologyView<MeshConfig> TriaCells<MeshConfig>::active_cell(
    const ActiveIdx idx) const
{
  return CellTopologyView<MeshConfig>(this, FlatIdx(m_active_to_linear_idx[idx.id()]));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::deactivate_nonleafs_and_update_active_cell_positions()
{
  Uint nb_active_cells = 0;

  for (Uint c = 0; c < m_status.size(); ++c)
  {
    if (m_cell_children[c].empty())
    {
      m_status[c] = EntityStatus::Active;
      nb_active_cells++;
    }
    else
    {
      m_status[c] = EntityStatus::NotActive;
    }
  }

  /*
  for (Uint c = 0; c < m_status.size(); ++c)
  {
    if (m_status[c] == EntityStatus::Active)
    {
      nb_active_cells++;
    }
    else
    {
      m_status[c] = EntityStatus::NotActive;
    }
  }
  */

  m_active_to_linear_idx.resize(nb_active_cells);
  nb_active_cells = 0;

  for (Uint c = 0; c < m_status.size(); ++c)
  {
    if (m_status[c] == EntityStatus::Active)
    {
      m_active_to_linear_idx[nb_active_cells] = c;
      nb_active_cells++;
    }
  }

  // Do not use 'resize' with initialization: this is going
  // to initialize the entries in m_linear_to_active_idx ONLY
  // if the new vector size is bigger than the previous one - and
  // even then, only the 'extra' entries will be initalized!
  // m_linear_to_active_idx.resize(m_status.size(), INVALID_CELL_ID);

  m_linear_to_active_idx.resize(m_status.size());
  m_linear_to_active_idx.assign(m_linear_to_active_idx.size(), INVALID_CELL_ID);

  for (Uint ac = 0; ac < m_active_to_linear_idx.size(); ++ac)
  {
    m_linear_to_active_idx[m_active_to_linear_idx[ac]] = ac;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> TriaCells<MeshConfig>::cell_parent(
    CellTopologyView<MeshConfig> const &cell) const
{
  const Uint parent_idx = m_cell_parent[cell.linear_pos_idx().id()];
  return CellTopologyView<MeshConfig>(this, FlatIdx(parent_idx));
}

// ----------------------------------------------------------------------------

/*
template <typename MeshConfig>
void CellConnectivity<MeshConfig>::set_parent(const Uint cell_idx, const Uint
parent_idx)
{
  m_cell_parent[cell_idx] = parent_idx;
}
*/

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<CellTopologyView<MeshConfig>> TriaCells<MeshConfig>::cell_children(
    const FlatIdx idx) const
{
  const std::vector<Uint> &child_ids = m_cell_children[idx.id()];

  std::vector<CellTopologyView<MeshConfig>> children(child_ids.size());

  for (Uint c = 0; c < child_ids.size(); ++c)
  {
    children[c] = CellTopologyView<MeshConfig>(this, FlatIdx(child_ids[c]));
  }

  // Hopefully move semantics applies here when returning a vector ...
  return children;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::path(const CellTopologyView<MeshConfig> &start_tcell,
                                 CellTopologyView<MeshConfig> &zero_level_cell,
                                 std::vector<Uint> &path_entries) const
{
  path_entries.resize(0);

  Uint start_tcell_lin_id = start_tcell.linear_pos_idx().id();

  if (m_level[start_tcell_lin_id] == 0)
  {
    zero_level_cell = start_tcell;
    return;
  }

  Uint zero_level_tcell_lin_id = start_tcell.linear_pos_idx().id();

  while (m_level[zero_level_tcell_lin_id] > 0)
  {
    const Uint parent_level_tcell_lin_id = m_cell_parent[zero_level_tcell_lin_id];
    const std::vector<Uint> &children    = m_cell_children[parent_level_tcell_lin_id];

    for (Uint loc_id = 0; loc_id < children.size(); ++loc_id)
    {
      if (children[loc_id] == zero_level_tcell_lin_id)
      {
        path_entries.push_back(loc_id);
        zero_level_tcell_lin_id = parent_level_tcell_lin_id;
        break;
      }
    }
  }

  zero_level_cell = cell(FlatIdx(zero_level_tcell_lin_id));
  std::reverse(path_entries.begin(), path_entries.end());
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const CellTopologyView<MeshConfig> TriaCells<MeshConfig>::path_leaf(const CellPath &path) const
{
  Uint curr_cell_lin_id = path.base_cell_id().id();
  for (Uint i = 0; i < path.size(); ++i)
  {
    const std::vector<Uint> &children = m_cell_children[curr_cell_lin_id];
    if (children.empty())
    {
      break;
    }
    curr_cell_lin_id = children[path.path_entry(i)];
  }

  return this->cell(FlatIdx(curr_cell_lin_id));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::add_child_indexes(const Uint cell_idx,
                                              std::vector<Uint> const &child_indexes,
                                              const EntityStatus &parent_status)
{
  m_status[cell_idx] = parent_status;
  // m_linear_to_active_idx[cell_idx] = INVALID_CELL_ID;
  for (Uint i = 0; i < child_indexes.size(); ++i)
  {
    m_cell_children[cell_idx].push_back(child_indexes[i]);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::mark_and_disconnect_children(const FlatIdx cell_idx,
                                                         const EntityStatus &parent_status,
                                                         const EntityStatus &child_status)
{
  const Uint cid = cell_idx.id();
  m_status[cid]  = parent_status;

  for (auto child_id : m_cell_children[cid])
  {
    m_status[child_id] = child_status;
  }

  m_cell_children[cid].clear();

  const ElemShape parent_shape = m_cell_type[cid].get().pt_set_id().elem_shape();
  m_cell_adapt_op[cid].change_type(parent_shape, CellTransform::NO_TRANS);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::remove_all_leaf_cells_w_status(const EntityStatus &status,
                                                           std::vector<FlatIdx> &new_cell_id)
{
  const Uint nb_old_cells = m_cell_type.size();

  if (nb_old_cells == 0)
  {
    std::cerr << "TriaCells: remove_all_leaf_cells_w_status: error, no cells "
                 "to remove"
              << std::endl;
    return;
  }

  // nb_of_precedent_cells_to_remove is a vector
  // which stores on position c the number of cells with
  // index STRICTLY SMALLER THAN c that will be removed
  new_cell_id.resize(nb_old_cells);
  new_cell_id.assign(nb_old_cells, FlatIdx(INVALID_CELL_ID));

  std::vector<bool> remove_cell;
  remove_cell.assign(nb_old_cells, false);
  Uint nb_cells_to_remove = 0;

  for (Uint c = 0; c < m_status.size(); ++c)
  {
    if ((m_status[c] == status) && (m_cell_children[c].empty()))
    {
      remove_cell[c] = true;
      nb_cells_to_remove++;
    }
    else
    {
      new_cell_id[c] = FlatIdx(c - nb_cells_to_remove);
    }
  }

  /*
  std::cout << "New cell ids: " << std::endl;
  for (auto id : new_cell_id)
  {
    std::cout << id.id() << " ";
  }
  std::cout << std::endl;

  std::cout << "Total number of cells to remove: " << nb_cells_to_remove <<
  std::endl;
  */

  std::vector<Uint> cell_ids_to_remove;
  cell_ids_to_remove.reserve(nb_cells_to_remove);
  cell_ids_to_remove.resize(0);

  // Go through a vector of cells to be removed
  // Store the ids of cells scheduled for removal in a vector
  // Remove all children of parents of cells scheduled for removal
  // Change the type of adaptation operation applied to parent to "DO_NOTHING'
  for (Uint c = 0; c < remove_cell.size(); ++c)
  {
    if (remove_cell[c])
    {
      cell_ids_to_remove.push_back(c);
      const Uint parent_id = m_cell_parent[c];
      m_cell_children[parent_id].resize(0);

      const ElemShape parent_shape = m_cell_type[parent_id].get().pt_set_id().elem_shape();
      m_cell_adapt_op[parent_id].change_type(parent_shape, CellTransform::NO_TRANS);
    }
  }

  // Position where data can be copied
  Uint fill_pos = 0;
  // Position to first entry that should be preserved
  Uint valid_data_pos = 0;

  Uint nb_new_cells = 0;

  // I) Modify all member data that are stored in blocks
  Uint nb_face_entries_after_removal = 0;
  for (Uint c = 0; c < nb_old_cells; ++c)
  {
    const Uint block_data_len = m_cell_face_offset[c + 1] - m_cell_face_offset[c];
    // If cell data should be preserved and there were
    // some data removed prior to processing the cell,
    // copy the data on the first available position
    if (!remove_cell[c])
    {
      nb_face_entries_after_removal += block_data_len;
      if (fill_pos != m_cell_face_offset[c])
      {
        for (Uint i = 0; i < block_data_len; ++i)
        {
          m_cell_faces[fill_pos] = m_cell_faces[valid_data_pos];
          fill_pos++;
          valid_data_pos++;
        }
      }
      else
      {
        fill_pos += block_data_len;
        valid_data_pos += block_data_len;
      }

      nb_new_cells++;
    }
    else
    {
      valid_data_pos += block_data_len;
    }
  }

  m_cell_faces.resize(nb_face_entries_after_removal);

  fill_pos = 0;

  // II) Modify all remaining information
  // Turn temporarily offsets into block sizes
  // Size of block c is stored on position [c+1]
  // m_cell_face_offset[0] will remain = 0
  for (Uint c = nb_old_cells; c > 0; --c)
  {
    m_cell_face_offset[c] = m_cell_face_offset[c] - m_cell_face_offset[c - 1];
  }

  m_cell_face_offset[0] = 0;

  // We work with the assumption that since cells are only being removed,
  // the inequality
  //              new_cell_id[c] <= c
  // is ALWAYS satisfied

  for (Uint c = 0; c < nb_old_cells; ++c)
  {
    if (!remove_cell[c])
    {

      std::vector<Uint> &old_cell_children = m_cell_children[c];

      if (fill_pos != c)
      {
        m_cell_type[fill_pos]     = m_cell_type[c];
        m_cell_adapt_op[fill_pos] = m_cell_adapt_op[c];
        m_level[fill_pos]         = m_level[c];
        m_status[fill_pos]        = m_status[c];

        // Shift the child vector
        std::vector<Uint> &new_cell_children = m_cell_children[fill_pos];

        new_cell_children.resize(old_cell_children.size());

        // Renumber the shifted child ids to take into account
        // the fact that some cells will be removed and so
        // the global cell numbering will be affected
        for (Uint i = 0; i < new_cell_children.size(); ++i)
        {
          new_cell_children[i] = new_cell_id[old_cell_children[i]].id();
        }
        old_cell_children.resize(0);

        m_cell_face_offset[fill_pos + 1] = m_cell_face_offset[c + 1];
      }
      else
      {
        // No need to shift the child cell vector, only
        // update the child cell ids
        for (Uint i = 0; i < old_cell_children.size(); ++i)
        {
          const Uint old_child_id = old_cell_children[i];
          old_cell_children[i]    = new_cell_id[old_child_id].id();
        }
      }

      fill_pos++;
    }

  } // Loop over old cells

  m_cell_type.resize(nb_new_cells);
  m_cell_adapt_op.resize(nb_new_cells);
  m_level.resize(nb_new_cells);
  m_status.resize(nb_new_cells);
  m_cell_children.resize(nb_new_cells);

  // Convert block sizes back to offsets
  for (Uint c = 1; c <= nb_new_cells; ++c)
  {
    m_cell_face_offset[c] += m_cell_face_offset[c - 1];
  }

  // Re-generate information about cell parents
  // For each cell that has children, inform all children
  // about the id of the parent

  // At the same time, make sure that leaf-level cells are active
  // and all other cells are marked as not active
  // Note that it is NOT sufficient to remove the data of cells
  // scheduled for removal and preserve all remaining information
  // Cells that have now become leaf-level cells DUE TO COARSENING
  // have to CHANGE their status from NotActive to Active

  m_cell_parent.resize(nb_new_cells);
  m_cell_parent.assign(nb_new_cells, INVALID_CELL_ID);

  for (Uint c = 0; c < m_cell_children.size(); ++c)
  {
    if (!m_cell_children[c].empty())
    {
      for (Uint i = 0; i < m_cell_children[c].size(); ++i)
      {
        const Uint child_id     = m_cell_children[c][i];
        m_cell_parent[child_id] = c;
      }
      m_status[c] = EntityStatus::NotActive;
    }
    else
    {
      m_status[c] = EntityStatus::Active;
    }
  }

  // Update indices referring to active cells and also storing
  // the flat indices of active cells
  m_active_to_linear_idx.resize(0);
  for (Uint c = 0; c < m_status.size(); ++c)
  {
    if (m_status[c] == EntityStatus::Active)
    {
      m_active_to_linear_idx.push_back(c);
    }
  }

  m_linear_to_active_idx.resize(m_status.size(), INVALID_CELL_ID);
  for (Uint ac = 0; ac < m_active_to_linear_idx.size(); ++ac)
  {
    m_linear_to_active_idx[m_active_to_linear_idx[ac]] = ac;
  }

  // Update coordinates array
  m_cell_coords.remove_blocks(cell_ids_to_remove);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::emplace_cells(std::vector<StdRegion> &cell_types,
                                          std::vector<Uint> &cell_faces,
                                          common::BlockArray<Real, Uint> const &cell_coordinates)
{
  m_max_level = 1;

  m_cell_type.swap(cell_types);
  m_cell_adapt_op.resize(m_cell_type.size());
  m_level.resize(m_cell_type.size(), 0);
  m_status.resize(m_cell_type.size(), EntityStatus::Active);
  m_active_to_linear_idx.resize(m_cell_type.size());
  for (Uint i = 0; i < m_cell_type.size(); ++i)
  {
    m_active_to_linear_idx[i] = i;
  }

  m_linear_to_active_idx.resize(m_status.size(), INVALID_CELL_ID);
  for (Uint ac = 0; ac < m_active_to_linear_idx.size(); ++ac)
  {
    m_linear_to_active_idx[m_active_to_linear_idx[ac]] = ac;
  }

  m_cell_faces.swap(cell_faces);

  for (Uint i = 0; i < m_cell_type.size(); ++i)
  {
    const PointSetTag original_tag = m_cell_type[i].get().pt_set_id();
    /*
    // const PointSetTag p1_tag(original_tag.elem_shape(), P1,
    original_tag.ref_topology()); const PointSetTag
    p1_tag(original_tag.elem_shape(), P1, PointSetID::Equidist);
    m_cell_type[i].change_type(p1_tag);
    */

    m_cell_adapt_op[i].change_type(original_tag.elem_shape(), CellTransform::NO_TRANS);
  }

  m_cell_face_offset.resize(m_cell_type.size() + 1);
  m_cell_face_offset[0] = 0;

  for (Uint i = 1; i < m_cell_face_offset.size(); ++i)
  {
    m_cell_face_offset[i] =
        m_cell_face_offset[i - 1] + m_cell_type[i - 1].get().nb_entities(MeshConfig::TDIM - 1);
  }

  m_cell_parent.resize(m_cell_type.size());
  std::fill(m_cell_parent.begin(), m_cell_parent.end(), INVALID_CELL_ID);
  m_cell_children.resize(m_cell_type.size());

  m_cell_coords = cell_coordinates;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename TraceTopologyIterator>
void TriaCells<MeshConfig>::rebuild_neighbour_information(const TraceTopologyIterator &facets_begin,
                                                          const TraceTopologyIterator &facets_end)
{
  // I) Generate facet offsets

  TraceTopologyIterator facet_it(facets_begin);

  // Number of facets of each cell
  // This will be later on turned into the vector of facet offsets
  std::vector<Uint> facet_offsets;
  facet_offsets.resize(m_cell_type.size() + 1);
  facet_offsets.assign(m_cell_type.size() + 1, 0);

  for (; facet_it != facets_end; ++facet_it)
  {
    const TraceIncidences facet_block = facet_it->incidences();
    if (all_incident_cells_have_status(facet_block, EntityStatus::Active))
    {
      for (Uint i = 0; i < facet_block.size(); ++i)
      {
        const Uint cell_idx = facet_block.cell_id(i);
        facet_offsets[cell_idx]++;
      }
    }
  }

  const Uint tot_nb_facets = std::accumulate(facet_offsets.begin(), facet_offsets.end(), 0);

  if (tot_nb_facets == 0)
  {
    std::cerr << "TriaCells::rebuild_neighbour_information: no facets." << std::endl;
    return;
  }

  // Now turn this array into the array of offsets.
  // First shift all data in the array to the right by one position
  // so that the first entry is 0
  for (Uint i = facet_offsets.size() - 1; i > 0; --i)
  {
    facet_offsets[i] = facet_offsets[i - 1];
  }

  facet_offsets[0] = 0;

  for (Uint i = 1; i < facet_offsets.size(); ++i)
  {
    facet_offsets[i] += facet_offsets[i - 1];
  }

  // II) Fill the facet data information

  std::vector<Uint> facet_ids;
  facet_ids.resize(tot_nb_facets);
  facet_ids.assign(tot_nb_facets, INVALID_CELL_ID);

  std::vector<SUint> local_ids;
  local_ids.resize(tot_nb_facets);

  for (facet_it = facets_begin; facet_it != facets_end; ++facet_it)
  {
    const TraceIncidences facet_block = facet_it->incidences();

    if (all_incident_cells_have_status(facet_block, EntityStatus::Active))
    {

      for (Uint i = 0; i < facet_block.size(); ++i)
      {
        const Uint cell_idx = facet_block.cell_id(i);
        for (Uint j = facet_offsets[cell_idx]; j < facet_offsets[cell_idx + 1]; ++j)
        {
          if (facet_ids[j] == INVALID_CELL_ID)
          {
            facet_ids[j] = facet_block.idx();
            local_ids[j] = facet_block.local_id(i);
            break;
          }
        }
      } // Loop over entries of one facet block
    }
  } // Iteration over facet blocks

  // III) Order facets of each cell so that they go according to
  //      increasing local id

  std::vector<Uint> one_cell_facets_reordered;

  for (Uint i = 0; i < (facet_offsets.size() - 1); ++i)
  {
    const Uint facet_block_len = facet_offsets[i + 1] - facet_offsets[i];

    if (facet_block_len > 0)
    {
      one_cell_facets_reordered.resize(0);

      const Uint pos_first = facet_offsets[i];
      const Uint pos_last  = facet_offsets[i + 1] - 1;

      SUint local_id_min = local_ids[pos_first];
      SUint local_id_max = local_ids[pos_first];

      for (Uint pos = pos_first + 1; pos <= pos_last; ++pos)
      {
        local_id_min = std::min(local_id_min, local_ids[pos]);
        local_id_max = std::max(local_id_max, local_ids[pos]);
      }

      for (SUint local_id = local_id_min; local_id <= local_id_max; ++local_id)
      {
        for (Uint pos = pos_first; pos <= pos_last; ++pos)
        {
          if (local_ids[pos] == local_id)
          {
            one_cell_facets_reordered.push_back(facet_ids[pos]);
          }
        }
      }

      // Copy the reordered facets back in the array 'facet_ids'
      std::copy(one_cell_facets_reordered.begin(), one_cell_facets_reordered.end(),
                facet_ids.data() + pos_first);

      /*
      for (Uint i = 0, pos = pos_begin; pos < pos_end; ++i, ++pos)
      {
        facet_ids[pos] = one_cell_facets_reordered[i];
      }
      */

    } // If the block length is > 0
  }   // Loop over all facet blocks

  std::swap(m_cell_face_offset, facet_offsets);
  std::swap(m_cell_faces, facet_ids);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename StdRegChangeRule>
void TriaCells<MeshConfig>::change_std_region_types(const StdRegChangeRule &rule)
{
  // ***********************************************************************
  // Compute the matrix which will enable interpolation to the new
  // nodes inside each element
  // ***********************************************************************

  sf::ShapeFunction shape_func;

  // This map associates a matrix of shape function values to an interpolation
  // point set id

  using interp_op_storage_t =
      std::map<std::tuple<PointSetTag, PointSetTag>, std::unique_ptr<math::DenseDMat<Real>>>;
  interp_op_storage_t type_to_interpolation_matrix;

  StdRegion new_std_reg;

  std::unique_ptr<std::vector<Real>> new_coords(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> new_offsets(new std::vector<Uint>());

  new_offsets->reserve(m_cell_type.size() + 1);
  new_offsets->resize(0);
  new_offsets->push_back(0);

  Uint nb_node_entries = 0;

  for (auto old_std_reg : m_cell_type)
  {
    const PointSetTag old_std_reg_tag = old_std_reg.get().pt_set_id();
    const PointSetTag new_std_reg_tag = rule(old_std_reg_tag);

    if (old_std_reg_tag.elem_shape() != new_std_reg_tag.elem_shape())
    {
      std::cerr << "TriaCells: can't upgrade a standard region by changing its "
                   "element shape!"
                << std::endl;
      return;
    }
    if (new_std_reg_tag.poly_order() < P1)
    {
      std::cerr << "TriaCells: can't upgrade tp a standard region with P = 0!" << std::endl;
      return;
    }

    new_std_reg.change_type(new_std_reg_tag);
    const Uint nb_new_interp_pts = new_std_reg.get().nb_nodes();

    const interp_op_storage_t::key_type key(old_std_reg_tag, new_std_reg_tag);
    interp_op_storage_t::iterator interp_it = type_to_interpolation_matrix.find(key);

    if (interp_it == type_to_interpolation_matrix.end())
    {
      const sf::SFTag old_sf_tag(old_std_reg_tag.elem_shape(), SFunc::Lagrange,
                                 old_std_reg_tag.poly_order(), ModalBasis::Modal);
      shape_func.change_type(old_std_reg_tag, old_sf_tag);

      StdRegion::value_type::coordinates_type const &new_std_reg_coords =
          new_std_reg.get().coordinates();

      const Uint nb_old_interp_pts = shape_func.get().nb_dof();

      std::unique_ptr<math::DenseDMat<Real>> interp_matrix(
          new math::DenseDMat<Real>(nb_new_interp_pts, nb_old_interp_pts));

      shape_func.get().compute_ref_values(new_std_reg_coords, *interp_matrix);
      // type_to_interpolation_matrix.insert(std::make_pair(key,
      // std::move(interp_matrix)));
      type_to_interpolation_matrix.emplace(key, std::move(interp_matrix));
    }

    nb_node_entries += MeshConfig::GDIM * nb_new_interp_pts;
    // Store the new block size in 'new_offsets'
    new_offsets->push_back(MeshConfig::GDIM * nb_new_interp_pts);
  }

  // Turn block sizes into offsets
  for (Uint i = 1; i < new_offsets->size(); ++i)
  {
    (*new_offsets)[i] += (*new_offsets)[i - 1];
  }

  new_coords->reserve(nb_node_entries);
  new_coords->resize(0);

  std::vector<Real> tmp_coord_old, tmp_coord_new;

  for (Uint i = 0; i < m_cell_type.size(); ++i)
  {
    const PointSetTag old_std_reg_tag = m_cell_type[i].get().pt_set_id();
    const PointSetTag new_std_reg_tag = rule(old_std_reg_tag);

    new_std_reg.change_type(new_std_reg_tag);

    const interp_op_storage_t::key_type key(old_std_reg_tag, new_std_reg_tag);
    interp_op_storage_t::const_iterator interp_it = type_to_interpolation_matrix.find(key);

    const math::DenseDMat<Real> &interp_op = *(interp_it->second);

    const Uint nb_nodes_old = interp_op.cols();
    const Uint nb_nodes_new = interp_op.rows();
    tmp_coord_old.resize(nb_nodes_old * MeshConfig::GDIM);
    tmp_coord_new.resize(nb_nodes_new * MeshConfig::GDIM);

    math::DenseMatView<Real> old_cell_coord(tmp_coord_old.data(), MeshConfig::GDIM, nb_nodes_old,
                                            MeshConfig::GDIM);
    math::DenseMatView<Real> new_cell_coord(tmp_coord_new.data(), MeshConfig::GDIM, nb_nodes_new,
                                            MeshConfig::GDIM);

    const common::ArrayView<const Real, _1D, Uint> coord_block = m_cell_coords.const_block(i);

    for (Uint n = 0; n < nb_nodes_old; ++n)
    {
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        old_cell_coord(n, d) = coord_block[n * MeshConfig::GDIM + d];
      }
    }

    new_cell_coord = interp_op * old_cell_coord;

    for (Uint n = 0; n < nb_nodes_new; ++n)
    {
      for (Uint d = 0; d < MeshConfig::GDIM; ++d)
      {
        new_coords->push_back(new_cell_coord(n, d));
      }
    }

    m_cell_type[i] = new_std_reg;
  }

  std::cout << "Number of old coordinate entries = " << m_cell_coords.size() << std::endl;
  std::cout << "Number of new coordinate entries = " << new_coords->size() << std::endl;
  m_cell_coords.build_from_offsets(std::move(new_coords), std::move(new_offsets));
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename GeoTransformRule>
void TriaCells<MeshConfig>::geo_transform(const GeoTransformRule &rule)
{
  const Uint GDIM = MeshConfig::GDIM;
  std::vector<Real> cell_coord_data_new;
  std::vector<Real> one_node_coord_data_old(GDIM);

  for (Uint c = 0; c < m_cell_coords.nb_blocks(); ++c)
  {
    // Copy existing cell coordinates into temporary storage
    common::ArrayView<const Real, _1D, Uint> cell_coord_old = m_cell_coords.const_block(c);

    const math::DenseConstVecView<Real> one_node_coord_old(one_node_coord_data_old.data(), GDIM);

    const Uint nb_elem_nodes = m_cell_type[c].get().nb_nodes();
    cell_coord_data_new.resize(nb_elem_nodes * GDIM);

    for (Uint n = 0; n < nb_elem_nodes; ++n)
    {
      math::DenseVecView<Real> one_node_coord_new(cell_coord_data_new.data() + n * GDIM, GDIM);
      for (Uint d = 0; d < GDIM; ++d)
      {
        one_node_coord_data_old[d] = cell_coord_old[n * GDIM + d];
      }
      rule(one_node_coord_old, one_node_coord_new);
      for (Uint d = 0; d < GDIM; ++d)
      {
        cell_coord_data_new[n * GDIM + d] = one_node_coord_new[d];
      }
    }
    const common::ArrayView<const Real, _1D, Uint> cell_coord_new(cell_coord_data_new.data(),
                                                                  cell_coord_data_new.size());
    m_cell_coords.insert_block(c, cell_coord_new);
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
template <typename GeoTransformRule>
void TriaCells<MeshConfig>::geo_transform(const ActiveIdx idx, const GeoTransformRule &rule)
{
  const Uint GDIM = MeshConfig::GDIM;
  std::vector<Real> cell_coord_data_new;
  std::vector<Real> one_node_coord_data_old(GDIM);

  const Uint flat_idx = m_active_to_linear_idx[idx.id()];

  // Copy existing cell coordinates into temporary storage
  common::ArrayView<const Real, _1D, Uint> cell_coord_old = m_cell_coords.const_block(flat_idx);

  const math::DenseConstVecView<Real> one_node_coord_old(one_node_coord_data_old.data(), GDIM);

  const Uint nb_elem_nodes = m_cell_type[flat_idx].get().nb_nodes();
  cell_coord_data_new.resize(nb_elem_nodes * GDIM);

  for (Uint n = 0; n < nb_elem_nodes; ++n)
  {
    math::DenseVecView<Real> one_node_coord_new(cell_coord_data_new.data() + n * GDIM, GDIM);
    for (Uint d = 0; d < GDIM; ++d)
    {
      one_node_coord_data_old[d] = cell_coord_old[n * GDIM + d];
    }

    rule(n, one_node_coord_old, one_node_coord_new);

    /*
    for (Uint d = 0; d < GDIM; ++d)
    {
      cell_coord_data_new[n * GDIM + d] = one_node_coord_new[d];
    }
    */
  }
  const common::ArrayView<const Real, _1D, Uint> cell_coord_new(cell_coord_data_new.data(),
                                                                cell_coord_data_new.size());

  m_cell_coords.insert_block(flat_idx, cell_coord_new);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::clear()
{
  m_max_level = 0;
  m_cell_faces.resize(0);

  m_cell_face_offset.resize(1);
  m_cell_face_offset[0] = 0;

  m_cell_type.resize(0);
  m_cell_adapt_op.resize(0);
  m_level.resize(0);
  m_status.resize(0);
  m_active_to_linear_idx.resize(0);
  m_linear_to_active_idx.resize(0);
  m_cell_parent.resize(0);
  m_cell_children.resize(0);
  m_cell_coords.clear();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::first_active_cell_pos() const
{
  return m_active_to_linear_idx[0];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::last_active_cell_pos() const
{
  return m_active_to_linear_idx.back();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
std::ostream &operator<<(std::ostream &os, const TriaCells<MeshConfig> &ml)
{
  os << "List of all cells in cell connectivity" << std::endl;
  for (Uint cell = 0; cell < ml.m_cell_type.size(); ++cell)
  {
    os << "[" << cell << "]: " << ml.cell_type(FlatIdx(cell)).get().pt_set_id().as_string()
       << std::endl;

    if (ml.m_level[cell] == 0u)
    {
      std::cout << "      Parent: None" << std::endl;
    }
    else
    {
      std::cout << "      Parent: " << ml.m_cell_parent[cell] << std::endl;
    }

    os << "      Children:" << std::endl;
    for (Uint c = 0; c < ml.m_cell_children[cell].size(); ++c)
    {
      os << "        " << ml.m_cell_children[cell][c] << std::endl;
    }
    os << "      Faces:" << std::endl;
    os << "        " << ml.cell_faces(FlatIdx(cell)) << std::endl;
    os << "      Face offsets: (" << ml.m_cell_face_offset[cell] << ","
       << ml.m_cell_face_offset[cell + 1] << ")" << std::endl;
    os << "      Status:" << ml.m_status[cell] << std::endl;
    os << "      Level:" << ml.m_level[cell] << std::endl;
    if (ml.m_status[cell] == EntityStatus::Active)
    {
      os << "      Active id:" << ml.m_linear_to_active_idx[cell] << std::endl;
    }
    os << "      Cell adapt operation = "
       << ml.cell_adapt_op(FlatIdx(cell)).get().cell_adapt_op_tag().as_string() << std::endl;
    os << std::endl;
  }

  os << "Linear to active index (total size " << ml.m_linear_to_active_idx.size()
     << ") : " << std::endl;
  for (auto idx : ml.m_linear_to_active_idx)
  {
    os << idx << " ";
  }
  os << std::endl;

  os << "Active to linear index (total size " << ml.m_active_to_linear_idx.size()
     << ") : " << std::endl;
  for (auto idx : ml.m_active_to_linear_idx)
  {
    os << idx << " ";
  }
  os << std::endl;

  return os;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::clone(const TriaCells &cell_conn_in, TriaCells &cell_conn_out)
{
  cell_conn_out.m_max_level = cell_conn_in.m_max_level;

  cell_conn_out.m_cell_faces.resize(cell_conn_in.m_cell_faces.size());
  std::copy(cell_conn_in.m_cell_faces.begin(), cell_conn_in.m_cell_faces.end(),
            cell_conn_out.m_cell_faces.begin());

  cell_conn_out.m_cell_face_offset.resize(cell_conn_in.m_cell_face_offset.size());
  std::copy(cell_conn_in.m_cell_face_offset.begin(), cell_conn_in.m_cell_face_offset.end(),
            cell_conn_out.m_cell_face_offset.begin());

  cell_conn_out.m_cell_type.resize(cell_conn_in.m_cell_type.size());
  std::copy(cell_conn_in.m_cell_type.begin(), cell_conn_in.m_cell_type.end(),
            cell_conn_out.m_cell_type.begin());

  cell_conn_out.m_cell_adapt_op.resize(cell_conn_in.m_cell_adapt_op.size());
  std::copy(cell_conn_in.m_cell_adapt_op.begin(), cell_conn_in.m_cell_adapt_op.end(),
            cell_conn_out.m_cell_adapt_op.begin());

  cell_conn_out.m_level.resize(cell_conn_in.m_level.size());
  std::copy(cell_conn_in.m_level.begin(), cell_conn_in.m_level.end(),
            cell_conn_out.m_level.begin());

  cell_conn_out.m_status.resize(cell_conn_in.m_status.size());
  std::copy(cell_conn_in.m_status.begin(), cell_conn_in.m_status.end(),
            cell_conn_out.m_status.begin());

  cell_conn_out.m_active_to_linear_idx.resize(cell_conn_in.m_active_to_linear_idx.size());
  std::copy(cell_conn_in.m_active_to_linear_idx.begin(), cell_conn_in.m_active_to_linear_idx.end(),
            cell_conn_out.m_active_to_linear_idx.begin());

  cell_conn_out.m_linear_to_active_idx.resize(cell_conn_in.m_linear_to_active_idx.size());
  std::copy(cell_conn_in.m_linear_to_active_idx.begin(), cell_conn_in.m_linear_to_active_idx.end(),
            cell_conn_out.m_linear_to_active_idx.begin());

  cell_conn_out.m_cell_parent.resize(cell_conn_in.m_cell_parent.size());
  std::copy(cell_conn_in.m_cell_parent.begin(), cell_conn_in.m_cell_parent.end(),
            cell_conn_out.m_cell_parent.begin());

  cell_conn_out.m_cell_children.resize(cell_conn_in.m_cell_children.size());

  for (Uint i = 0; i < cell_conn_in.m_cell_children.size(); ++i)
  {
    cell_conn_out.m_cell_children[i].resize(cell_conn_in.m_cell_children[i].size());
    std::copy(cell_conn_in.m_cell_children[i].begin(), cell_conn_in.m_cell_children[i].end(),
              cell_conn_out.m_cell_children[i].begin());
  }

  cell_conn_in.m_cell_coords = cell_conn_out.m_cell_coords;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::print_cell_coords_to_file(const std::string &filename) const
{
  std::ofstream outfile;
  outfile.open(filename.c_str());

  outfile.precision(10);

  for (Uint b = 0; b < m_cell_coords.nb_blocks(); ++b)
  {
    const common::ArrayView<const Real, _1D, Uint> cell_coords = m_cell_coords.const_block(b);
    const Uint nb_nodes = cell_coords.size() / MeshConfig::GDIM;

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint c = 0; c < MeshConfig::GDIM; ++c)
      {
        outfile << cell_coords[n * MeshConfig::GDIM + c] << " ";
      }
      outfile << std::endl;
    }
  }

  outfile.close();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const common::ArrayView<const Uint, _1D, Uint> TriaCells<MeshConfig>::cell_faces(
    const FlatIdx cell_idx) const
{
  const auto pos = cell_idx.id();
  const common::ArrayView<const Uint, _1D, Uint> faces(
      m_cell_faces.data() + m_cell_face_offset[pos],
      m_cell_face_offset[pos + 1] - m_cell_face_offset[pos]);
  return faces;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::reorder_cell_facets(const CellTopologyView<MeshConfig> &cell,
                                                const std::vector<Uint> &reordering)
{
  const Uint cell_idx = cell.linear_pos_idx().id();
  common::ArrayView<Uint, _1D, Uint> faces(m_cell_faces.data() + m_cell_face_offset[cell_idx],
                                           m_cell_face_offset[cell_idx + 1] -
                                               m_cell_face_offset[cell_idx]);

  std::vector<Uint> new_faces(faces.size());
  for (Uint i = 0; i < faces.size(); ++i)
  {
    new_faces[reordering[i]] = faces[i];
  }

  for (Uint i = 0; i < faces.size(); ++i)
  {
    faces[i] = new_faces[i];
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
StdRegion TriaCells<MeshConfig>::cell_type(const FlatIdx cell_idx) const
{
  return m_cell_type[cell_idx.id()];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const adapt::CellAdaptOp TriaCells<MeshConfig>::cell_adapt_op(const FlatIdx cell_idx) const
{
  return m_cell_adapt_op[cell_idx.id()];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const ActiveIdx TriaCells<MeshConfig>::active_idx(const FlatIdx idx) const
{
  return ActiveIdx(m_linear_to_active_idx[idx.id()]);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::nb_children(const FlatIdx idx) const
{
  return m_cell_children[idx.id()].size();
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
inline const mesh::CellGeometry<MeshConfig::GDIM> TriaCells<MeshConfig>::cell_geometry(
    const FlatIdx idx) const
{
  const std::tuple<const Real *, const Uint> raw_blk_dta = m_cell_coords.raw_data_block(idx.id());
  const math::DenseVecView<const Real> coord_view(std::get<0>(raw_blk_dta),
                                                  std::get<1>(raw_blk_dta));
  const mesh::CellGeometry<MeshConfig::GDIM> cell_geo(
      coord_view, m_cell_type[idx.id()].get().elem_entity(MeshConfig::TDIM, 0));
  return cell_geo;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const mesh::CellGeometry<MeshConfig::GDIM> TriaCells<MeshConfig>::cell_geometry(
    const FlatIdx idx, const Uint sub_entity_dim, const Uint local_id) const
{
  const std::tuple<const Real *, const Uint> raw_blk_dta = m_cell_coords.raw_data_block(idx.id());
  const math::DenseVecView<const Real> coord_view(std::get<0>(raw_blk_dta),
                                                  std::get<1>(raw_blk_dta));
  const std::shared_ptr<StdRegionEntity const> sub_entity =
      m_cell_type[idx.id()].get().elem_entity(sub_entity_dim, local_id);

  const mesh::CellGeometry<MeshConfig::GDIM> cell_geo(coord_view, sub_entity);
  return cell_geo;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::set_cell_adapt_op(const FlatIdx cell_idx,
                                              const adapt::CellAdaptOp &adapt_op)
{
  m_cell_adapt_op[cell_idx.id()] = adapt_op;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
Uint TriaCells<MeshConfig>::refinement_level(const FlatIdx cell_idx) const
{
  return m_level[cell_idx.id()];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
EntityStatus TriaCells<MeshConfig>::cell_status(const FlatIdx cell_idx) const
{
  return m_status[cell_idx.id()];
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
const std::vector<EntityStatus> &TriaCells<MeshConfig>::cell_status() const
{
  return m_status;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void TriaCells<MeshConfig>::change_status(const FlatIdx linear_pos_idx, const EntityStatus status)
{
  m_status[linear_pos_idx.id()] = status;
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
bool TriaCells<MeshConfig>::all_incident_cells_have_status(const TraceIncidences &incidences,
                                                           const EntityStatus status) const
{
  for (Uint i = 0; i < incidences.size(); ++i)
  {
    const Uint cell_id = incidences.cell_id(i);
    if (m_status[cell_id] != status)
    {
      return false;
    }
  }
  return true;
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace mesh

} // namespace pdekit

#endif
