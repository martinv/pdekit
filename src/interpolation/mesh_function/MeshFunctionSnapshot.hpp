#ifndef PDEKIT_Interpolation_Mesh_Function_Snapshot_hpp
#define PDEKIT_Interpolation_Mesh_Function_Snapshot_hpp

#include <map>
#include <vector>

#include "interpolation/ElementAdaptInterpolator.hpp"
#include "interpolation/mesh_function/ScalarMeshFunctionBase.hpp"
#include "interpolation/mesh_function/VectorMeshFunctionBase.hpp"
#include "math/DenseMatView.hpp"
#include "mesh/Tria.hpp"
#include "mesh/containers/CellPath.hpp"
#include "mesh/containers/DofMap.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

namespace detail
{

template <typename T>
class LocalElemPInterpolation
{
  public:
  /// Constructor
  LocalElemPInterpolation();

  /// Deleted copy constructor
  LocalElemPInterpolation(const LocalElemPInterpolation &other) = delete;

  /// ~Destructor
  ~LocalElemPInterpolation();

  /// Deleted assignment operator
  LocalElemPInterpolation &operator=(const LocalElemPInterpolation &rhs) = delete;

  /// Reset all internal data
  void reset();

  /// Compute interpolation
  void interpolate(const mesh::PointSetTag tag_elem_in, const mesh::PointSetTag tag_elem_out,
                   const math::DenseDMat<T> &values_in);

  /// Get interpolated values
  math::DenseConstMatView<T> interpolated_values() const;

  private:
  typedef std::pair<mesh::PointSetTag, mesh::PointSetTag> map_key_type;
  typedef std::unique_ptr<math::DenseDMat<T>> map_value_type;
  typedef std::map<map_key_type, map_value_type> operator_map_type;

  // Add new interpolation operator
  typename operator_map_type::iterator add_operator(map_key_type const &key);

  operator_map_type m_interp_operators;

  // Matrix representing the shape of interpolated data
  math::DenseMatView<T> m_result_proxy;

  // Raw interpolated data stored in array
  std::vector<T> m_work_data;
};

template <typename T>
LocalElemPInterpolation<T>::LocalElemPInterpolation()
{
}

template <typename T>
LocalElemPInterpolation<T>::~LocalElemPInterpolation()
{
}

template <typename T>
void LocalElemPInterpolation<T>::reset()
{
  m_interp_operators.clear();
}

template <typename T>
void LocalElemPInterpolation<T>::interpolate(const mesh::PointSetTag tag_elem_in,
                                             const mesh::PointSetTag tag_elem_out,
                                             const math::DenseDMat<T> &values_in)
{
  const map_key_type search_key(tag_elem_in, tag_elem_out);
  typename operator_map_type::iterator it = m_interp_operators.find(search_key);

  if (it == m_interp_operators.end())
  {
    it = add_operator(search_key);
  }

  const math::DenseDMat<T> &interp_op = (*it->second);

  const Uint nb_rows_result = interp_op.rows();
  const Uint nb_cols_result = values_in.cols();

  m_work_data.resize(nb_rows_result * nb_cols_result);
  m_result_proxy =
      math::DenseMatView<Real>(m_work_data.data(), nb_cols_result, nb_rows_result, nb_cols_result);
  m_result_proxy = interp_op * values_in;
}

template <typename T>
typename LocalElemPInterpolation<T>::operator_map_type::iterator LocalElemPInterpolation<
    T>::add_operator(map_key_type const &key)
{
  mesh::StdRegion elem_out;
  elem_out.change_type(key.second);

  const mesh::sf::SFTag sf_tag(key.first.elem_shape(), SFunc::Lagrange, key.first.poly_order(),
                               ModalBasis::Modal);

  mesh::sf::ShapeFunction sf;
  sf.change_type(key.first, sf_tag);

  map_value_type op_values = std::unique_ptr<math::DenseDMat<T>>(new math::DenseDMat<T>());

  op_values->resize(elem_out.get().nb_nodes(), sf.get().nb_dof());

  sf.get().compute_ref_values(elem_out.get().coordinates(), *op_values);

  const std::pair<typename operator_map_type::iterator, bool> res =
      m_interp_operators.insert(std::make_pair(key, std::move(op_values)));

  return res.first;
}

template <typename T>
math::DenseConstMatView<T> LocalElemPInterpolation<T>::interpolated_values() const
{
  const math::DenseConstMatView<T> values(m_work_data.data(), m_result_proxy.cols(),
                                          m_result_proxy.rows(), m_result_proxy.cols());
  return values;
}

} // namespace detail

// ----------------------------------------------------------------------------

template <typename T>
class MeshFunctionSnapshot
{
  public:
  /// Constructor
  MeshFunctionSnapshot();

  /// Deleted copy constructor
  MeshFunctionSnapshot(const MeshFunctionSnapshot &other) = delete;

  /// Destructor
  ~MeshFunctionSnapshot();

  /// Deleted assignment operator
  MeshFunctionSnapshot &operator=(const MeshFunctionSnapshot &rhs) = delete;

  /// Clear all internal data
  void clear();

  /// Create a snapshot of scalar function
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void create(const typename mesh::Tria<MeshConfig> &mesh_topology,
              const typename result_of::dof_map_t<MeshConfig> &dof_handler,
              const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  /// Create a snapshot of scalar function
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void create(const typename mesh::Tria<MeshConfig> &mesh_topology,
              const typename result_of::dof_map_t<MeshConfig> &dof_handler,
              const interpolation::VectorMeshFunctionBase<MFT> &mesh_function);

  /// Create a snapshot of scalar function which stores one value per active
  /// cell
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void create_cellwise(const typename mesh::Tria<MeshConfig> &mesh_topology,
                       const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                       const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  /// Restore a scalar function from snapshot
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void restore_function(const typename mesh::Tria<MeshConfig> &mesh_topology,
                        const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                        interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  /// Restore a scalar function from snapshot
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void restore_function_cellwise(const typename mesh::Tria<MeshConfig> &mesh_topology,
                                 const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                                 interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  /// Restore a scalar function from snapshot
  /// @param MFT ... Mesh Function Type
  template <typename MeshConfig, typename MFT>
  void restore_function(const typename mesh::Tria<MeshConfig> &mesh_topology,
                        const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                        interpolation::VectorMeshFunctionBase<MFT> &mesh_function);

  private:
  /// TYPES
  struct SnapshotData
  {
    /// Values of stored mesh function
    /// For each element, the values are stored by concatenating all values
    /// in all DOFs of the element corresponding to the first field,
    /// then all vertex values corresponding to the second field,
    /// then all vertex valus corresponding to the third field and so on
    ///
    /// E_i-F_j ... vector of field values in all vertices of element i,
    /// j-th field component Note that E_i-N_j is actually series of field
    /// values that has as entries as there are vertices in element!
    ///
    /// m_values has therefore the following structure
    /// [{ E_0-F_0, E_0-F_1, E_0-F_2 ..., E_0-F_m }, ..., { E_n-F_0,
    /// E_n-F_1,
    /// ..., E_n-F_m }]
    std::vector<T> values;

    /// Offsets delimiting mesh function values
    /// belonging to different elements/blocks
    std::vector<Uint> offsets;

    /// Linear cell index of each cell on which
    /// the function was defined
    std::vector<mesh::CellPath> cell_paths;

    /// Type of each element
    std::vector<mesh::StdRegion> cell_type;

    /// Number of fields in the function
    Uint nb_fields;
  };

  /// METHODS
  template <typename MeshConfig, typename MFT>
  void create_scalar_snapshot(const typename mesh::Tria<MeshConfig> &mesh_topology,
                              const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                              const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function,
                              SnapshotData &snapshot);

  template <typename MeshConfig, typename MFT>
  void create_cellwise_scalar_snapshot(
      const typename mesh::Tria<MeshConfig> &mesh_topology,
      const typename result_of::dof_map_t<MeshConfig> &dof_handler,
      const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function, SnapshotData &snapshot);

  template <typename MeshConfig, typename MFT>
  void create_vector_snapshot(const typename mesh::Tria<MeshConfig> &mesh_topology,
                              const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                              const interpolation::VectorMeshFunctionBase<MFT> &mesh_function,
                              SnapshotData &snapshot);

  template <typename MeshConfig, typename MFT>
  void restore_scalar_function(const typename mesh::Tria<MeshConfig> &mesh_topology,
                               const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                               const SnapshotData &snapshot,
                               interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  template <typename MeshConfig, typename MFT>
  void restore_scalar_function_cellwise(
      const typename mesh::Tria<MeshConfig> &input_mesh,
      const typename result_of::dof_map_t<MeshConfig> &dof_handler, const SnapshotData &snapshot,
      interpolation::ScalarMeshFunctionBase<MFT> &mesh_function);

  template <typename MeshConfig, typename MFT>
  void restore_vector_function(const typename mesh::Tria<MeshConfig> &mesh_topology,
                               const typename result_of::dof_map_t<MeshConfig> &dof_handler,
                               const SnapshotData &snapshot,
                               interpolation::VectorMeshFunctionBase<MFT> &mesh_function);

  // typedef std::map<std::string, std::unique_ptr<SnapshotData>>
  // snapshot_storage;
  using snapshot_storage = std::map<std::string, std::unique_ptr<SnapshotData>>;

  snapshot_storage m_snapshots;
};

// ----------------------------------------------------------------------------

template <typename T>
MeshFunctionSnapshot<T>::MeshFunctionSnapshot()
{
}

// ----------------------------------------------------------------------------

template <typename T>
MeshFunctionSnapshot<T>::~MeshFunctionSnapshot()
{
  clear();
}

// ----------------------------------------------------------------------------

template <typename T>
void MeshFunctionSnapshot<T>::clear()
{
  m_snapshots.clear();
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it != m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::create (scalar function snapshot): " << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name() << " already exists."
              << std::endl;
    return;
  }

  std::unique_ptr<SnapshotData> snapshot_ptr(new SnapshotData);
  create_scalar_snapshot<MeshConfig, MFT>(mesh_topology, dof_handler, mesh_function, *snapshot_ptr);

  m_snapshots.insert(std::pair<std::string, std::unique_ptr<SnapshotData>>(
      wrapped_function.name(), std::move(snapshot_ptr)));
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it != m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::create (vector function snapshot): " << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name() << " already exists."
              << std::endl;
    return;
  }

  std::unique_ptr<SnapshotData> snapshot_ptr(new SnapshotData);
  create_vector_snapshot<MeshConfig, MFT>(mesh_topology, dof_handler, mesh_function, *snapshot_ptr);

  m_snapshots.insert(std::pair<std::string, std::unique_ptr<SnapshotData>>(
      wrapped_function.name(), std::move(snapshot_ptr)));
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create_cellwise(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it != m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::create (scalar function snapshot): " << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name() << " already exists."
              << std::endl;
    return;
  }

  std::unique_ptr<SnapshotData> snapshot_ptr(new SnapshotData);
  create_cellwise_scalar_snapshot<MeshConfig, MFT>(mesh_topology, dof_handler, mesh_function,
                                                   *snapshot_ptr);

  m_snapshots.insert(std::pair<std::string, std::unique_ptr<SnapshotData>>(
      wrapped_function.name(), std::move(snapshot_ptr)));
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_function(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it == m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::restore_function (scalar function "
                 "snapshot): "
              << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name()
              << " doesn't exist. Can not restore function." << std::endl;
    return;
  }

  restore_scalar_function<MeshConfig, MFT>(mesh_topology, dof_handler, *(it->second),
                                           mesh_function);
  m_snapshots.erase(it);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_function_cellwise(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it == m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::restore_function (scalar function "
                 "snapshot): "
              << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name()
              << " doesn't exist. Can not restore function." << std::endl;
    return;
  }

  restore_scalar_function_cellwise<MeshConfig, MFT>(mesh_topology, dof_handler, *(it->second),
                                                    mesh_function);
  m_snapshots.erase(it);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_function(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    interpolation::VectorMeshFunctionBase<MFT> &mesh_function)
{
  MFT const &wrapped_function            = mesh_function.wrapped_type();
  typename snapshot_storage::iterator it = m_snapshots.find(wrapped_function.name());

  if (it == m_snapshots.end())
  {
    std::cerr << "MeshFunctionSnapshot::restore_function (vector function "
                 "snapshot): "
              << std::endl;
    std::cerr << "Snapshot with name " << wrapped_function.name()
              << " doesn't exist. Can not restore function." << std::endl;
    return;
  }

  restore_vector_function<MeshConfig, MFT>(mesh_topology, dof_handler, *(it->second),
                                           mesh_function);
  m_snapshots.erase(it);
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create_scalar_snapshot(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function, SnapshotData &snapshot)
{
  snapshot.offsets.resize(dof_handler.nb_active_cells() + 1);
  snapshot.cell_paths.resize(dof_handler.nb_active_cells());
  snapshot.cell_type.resize(dof_handler.nb_active_cells());
  snapshot.nb_fields = 1;

  mesh::StdRegion std_region;

  // Number of entries to save
  Uint nb_entries = 0;

  mesh::CellTopologyView<MeshConfig> root_tcell;
  std::vector<Uint> path_entries;

  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(c));
    snapshot.offsets[c + 1]     = cell.nb_vert();
    nb_entries += cell.nb_vert();

    const mesh::CellTopologyView<MeshConfig> tcell = mesh_topology.active_cell(mesh::ActiveIdx(c));

    mesh_topology.path(tcell, root_tcell, path_entries);
    snapshot.cell_paths[c].build(mesh::FlatIdx(root_tcell.linear_pos_idx()), path_entries);

    std_region.change_type(cell.pt_set_id());
    snapshot.cell_type[c] = std_region;
  }

  snapshot.values.resize(nb_entries);

  for (Uint i = 1; i < snapshot.offsets.size(); ++i)
  {
    snapshot.offsets[i] += snapshot.offsets[i - 1];
  }

  // Finally, cache the function values
  MFT const &function_wrapped = mesh_function.wrapped_type();

  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(c));

    Uint init_pos = snapshot.offsets[c];
    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      snapshot.values[init_pos++] = function_wrapped[cell.vertex(v)];
    }
  }
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create_cellwise_scalar_snapshot(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::ScalarMeshFunctionBase<MFT> &mesh_function, SnapshotData &snapshot)
{
  using tcell_t = mesh::CellTopologyView<MeshConfig>;

  snapshot.values.reserve(dof_handler.nb_active_cells());
  snapshot.offsets.reserve(dof_handler.nb_active_cells() + 1);
  snapshot.cell_paths.reserve(dof_handler.nb_active_cells());
  snapshot.cell_type.reserve(dof_handler.nb_active_cells());

  snapshot.values.resize(0);
  snapshot.offsets.resize(0);
  snapshot.offsets.push_back(0);
  snapshot.cell_paths.resize(0);
  snapshot.cell_type.resize(0);

  snapshot.nb_fields = 1;

  mesh::StdRegion std_region;

  mesh::CellTopologyView<MeshConfig> root_tcell;
  std::vector<Uint> path_entries;
  mesh::CellPath cpath;

  std::vector<bool> cell_processed(mesh_topology.nb_all_cells_in_all_levels(), false);

  MFT const &function_wrapped = mesh_function.wrapped_type();

  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const tcell_t tcell = mesh_topology.active_cell(mesh::ActiveIdx(c));

    if (!cell_processed[tcell.linear_pos_idx().id()])
    {
      T store_value{};

      // Check whether the cell went through 'green', i.e. anisotropic
      // refinement
      if (tcell.refinement_level() > 0)
      {
        const tcell_t parent_tcell = tcell.parent();
        const mesh::CellTransform parent_trans =
            parent_tcell.cell_adapt_op().get().cell_adapt_op_tag().adapt_op_id();

        // If the cell is green, then we cache the parent cell instead.
        // The function value cached is the 'average' value of all
        // children of the parent
        if (mesh::CellTransformTraits::is_aniso_refinement(parent_trans))
        {
          const std::vector<tcell_t> children_tcells = parent_tcell.children();
          store_value = function_wrapped[children_tcells[0].active_idx().id()];

          cell_processed[children_tcells[0].linear_pos_idx().id()] = true;

          for (Uint i = 1; i < children_tcells.size(); ++i)
          {
            const tcell_t &child                        = children_tcells[i];
            const T new_value                           = function_wrapped[child.active_idx().id()];
            store_value                                 = std::max(store_value, new_value);
            cell_processed[child.linear_pos_idx().id()] = true;
          }

          mesh_topology.path(parent_tcell, root_tcell, path_entries);
          cpath.build(mesh::FlatIdx(root_tcell.linear_pos_idx()), path_entries);

          const mesh::MeshEntity cell =
              dof_handler.cell(mesh::FlatIdx(parent_tcell.linear_pos_idx()));
          std_region.change_type(cell.pt_set_id());
        } // If this cell was refined anisotropically
        else
        {
          store_value = function_wrapped[tcell.active_idx().id()];
          mesh_topology.path(tcell, root_tcell, path_entries);
          cpath.build(mesh::FlatIdx(root_tcell.linear_pos_idx()), path_entries);

          const mesh::MeshEntity cell = dof_handler.cell(tcell.linear_pos_idx());
          std_region.change_type(cell.pt_set_id());

          cell_processed[tcell.linear_pos_idx().id()] = true;
        }

      } // If this is cell with refinement level > 0
      else
      {
        store_value = function_wrapped[tcell.active_idx().id()];
        cpath.build(mesh::FlatIdx(tcell.linear_pos_idx()), {});

        const mesh::MeshEntity cell = dof_handler.cell(tcell.linear_pos_idx());
        std_region.change_type(cell.pt_set_id());
        cell_processed[tcell.linear_pos_idx().id()] = true;
      }

      snapshot.values.push_back(store_value);
      snapshot.offsets.push_back(1);
      snapshot.cell_paths.push_back(cpath);
      snapshot.cell_type.push_back(std_region);

    } // If cell not yet processed

  } // Loop over active cells

  for (Uint i = 1; i < snapshot.offsets.size(); ++i)
  {
    snapshot.offsets[i] += snapshot.offsets[i - 1];
  }
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::create_vector_snapshot(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler,
    const interpolation::VectorMeshFunctionBase<MFT> &mesh_function, SnapshotData &snapshot)
{
  snapshot.offsets.resize(dof_handler.nb_active_cells() + 1);
  snapshot.cell_paths.resize(dof_handler.nb_active_cells());
  snapshot.cell_type.resize(dof_handler.nb_active_cells());

  MFT const &function_wrapped = mesh_function.wrapped_type();

  snapshot.nb_fields = function_wrapped.nb_fields();

  mesh::StdRegion std_region;

  // Number of entries to save
  Uint nb_entries = 0;

  mesh::CellTopologyView<MeshConfig> root_tcell;
  std::vector<Uint> path_entries;

  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(c));
    snapshot.offsets[c + 1]     = cell.nb_vert() * snapshot.nb_fields;
    nb_entries += cell.nb_vert() * snapshot.nb_fields;

    const mesh::CellTopologyView<MeshConfig> tcell = mesh_topology.active_cell(mesh::ActiveIdx(c));

    mesh_topology.path(tcell, root_tcell, path_entries);
    snapshot.cell_paths[c].build(mesh::FlatIdx(root_tcell.linear_pos_idx()), path_entries);

    std_region.change_type(cell.pt_set_id());
    snapshot.cell_type[c] = std_region;
  }

  snapshot.values.resize(nb_entries);

  for (Uint i = 1; i < snapshot.offsets.size(); ++i)
  {
    snapshot.offsets[i] += snapshot.offsets[i - 1];
  }

  // Finally, cache the function values
  for (Uint c = 0; c < dof_handler.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(c));

    const Uint init_pos        = snapshot.offsets[c];
    const Uint nb_vert_in_cell = cell.nb_vert();
    for (Uint v = 0; v < nb_vert_in_cell; ++v)
    {
      typename MFT::const_entry_type const node_value =
          function_wrapped.const_value(cell.vertex(v));

      for (Uint f = 0; f < snapshot.nb_fields; ++f)
      {
        snapshot.values[init_pos + nb_vert_in_cell * f + v] = node_value[f];
      }
    }
  }
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_scalar_function(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler, const SnapshotData &snapshot,
    interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  mesh::StdRegion std_region;
  Uint nb_entries = 0;

  for (Uint i = 0; i < snapshot.cell_paths.size(); ++i)
  {
    const mesh::CellTopologyView<MeshConfig> tcell =
        mesh_topology.path_leaf(snapshot.cell_paths[i]);
    if (tcell.status() == mesh::EntityStatus::Active)
    {
      const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(tcell.active_idx()));
      nb_entries += cell.nb_vert();
      // nb_entries += snapshot.offsets[i + 1] - snapshot.offsets[i];
    }
    else
    {
      const std::vector<mesh::CellTopologyView<MeshConfig>> cell_children = tcell.children();

      for (Uint c = 0; c < cell_children.size(); ++c)
      {
        const mesh::MeshEntity child_cell =
            dof_handler.active_cell(mesh::ActiveIdx(cell_children[c].active_idx()));
        nb_entries += child_cell.nb_vert();
      }
    }
  }

  MFT &function_wrapped = mesh_function.wrapped_type();
  function_wrapped.resize(nb_entries);

  ElementAdaptInterpolator local_h_interpolator;
  detail::LocalElemPInterpolation<T> local_p_interpolator;

  math::DenseDMat<T> parent_values;

  for (Uint i = 0; i < snapshot.cell_paths.size(); ++i)
  {
    const mesh::CellTopologyView<MeshConfig> tcell =
        mesh_topology.path_leaf(snapshot.cell_paths[i]);

    // CASE I: if the element remains active, then it was not adapted, or it
    // was p-adapted
    if (tcell.status() == mesh::EntityStatus::Active)
    {
      const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(tcell.active_idx()));

      const mesh::PointSetTag old_cell_tag = snapshot.cell_type[i].get().pt_set_id();

      // If the tag of the cell is unchanged, then the cell was not
      // adapted (in any way)
      if (old_cell_tag == cell.pt_set_id())
      {
        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          function_wrapped[cell.vertex(v)] = snapshot.values[snapshot.offsets[i] + v];
        }
      }
      // Else it was p-adapted
      else
      {
        // Put the cached values (before adaptation) into a matrix
        const Uint nb_parent_nodes = snapshot.cell_type[i].get().nb_nodes();
        parent_values.resize(nb_parent_nodes, snapshot.nb_fields);

        Uint cached_value_idx = snapshot.offsets[i];

        for (Uint n = 0; n < nb_parent_nodes; ++n)
        {
          parent_values(n, 0) = snapshot.values[cached_value_idx++];
        }

        local_p_interpolator.interpolate(old_cell_tag, cell.pt_set_id(), parent_values);
        const math::DenseConstMatView<T> interpolated_values =
            local_p_interpolator.interpolated_values();

        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          function_wrapped[cell.vertex(v)] = interpolated_values(v, 0);
        }
      }
    }
    // CASE II: if the element is no longer active, then it was h-adapted
    else
    {
      // Put the cached values (before adaptation) into a matrix
      const Uint nb_parent_nodes = snapshot.cell_type[i].get().nb_nodes();
      parent_values.resize(nb_parent_nodes, snapshot.nb_fields);

      Uint cached_value_idx = snapshot.offsets[i];

      for (Uint n = 0; n < nb_parent_nodes; ++n)
      {
        parent_values(n, 0) = snapshot.values[cached_value_idx++];
      }

      // Interpolate the (old) parent values on children
      local_h_interpolator.compute_child_values(snapshot.cell_type[i].get().pt_set_id(),
                                                tcell.cell_adapt_op().get().cell_adapt_op_tag(),
                                                parent_values);

      // Copy the values into the function
      const std::vector<mesh::CellTopologyView<MeshConfig>> cell_children = tcell.children();

      for (Uint c = 0; c < cell_children.size(); ++c)
      {
        const mesh::MeshEntity child_cell =
            dof_handler.active_cell(mesh::ActiveIdx(cell_children[c].active_idx()));
        const math::DenseConstMatView<T> child_values =
            local_h_interpolator.interpolated_child_values(c);
        for (Uint v = 0; v < child_cell.nb_vert(); ++v)
        {
          function_wrapped[child_cell.vertex(v)] = child_values(v, 0);
        }
      } // Loop over cell children
    }   // else
  }     // loop over cell_paths
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_scalar_function_cellwise(
    const typename mesh::Tria<MeshConfig> &input_mesh,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler, const SnapshotData &snapshot,
    interpolation::ScalarMeshFunctionBase<MFT> &mesh_function)
{
  using tcell_t = mesh::CellTopologyView<MeshConfig>;

  const Uint nb_active_cells = input_mesh.nb_active_cells();
  MFT &function_wrapped      = mesh_function.wrapped_type();

  function_wrapped.resize(nb_active_cells);
  function_wrapped.fill(T());

  std::queue<mesh::FlatIdx> elem_queue;

  for (Uint i = 0; i < snapshot.cell_paths.size(); ++i)
  {
    const T restore_value = snapshot.values[i];

    const tcell_t tcell = input_mesh.path_leaf(snapshot.cell_paths[i]);
    if (tcell.status() == mesh::EntityStatus::Active)
    {
      function_wrapped[tcell.active_idx().id()] = restore_value;
    }
    else
    {
      const std::vector<tcell_t> children = tcell.children();
      for (const tcell_t &child : children)
      {
        elem_queue.push(mesh::FlatIdx(child.linear_pos_idx()));
      }

      while (!elem_queue.empty())
      {
        const mesh::FlatIdx cid = elem_queue.front();
        elem_queue.pop();

        const tcell_t tmp_cell = input_mesh.cell(mesh::FlatIdx(cid.id()));

        if (tmp_cell.status() == mesh::EntityStatus::Active)
        {
          function_wrapped[tmp_cell.active_idx().id()] = restore_value;
        }
        else
        {
          const std::vector<tcell_t> tmp_children = tmp_cell.children();
          for (const tcell_t &tmp_child : tmp_children)
          {
            elem_queue.push(mesh::FlatIdx(tmp_child.linear_pos_idx()));
          }
        }
      } // while queue not empty

    } // if tcell is not active
  }   // loop over snapshot data
}

// ----------------------------------------------------------------------------

template <typename T>
template <typename MeshConfig, typename MFT>
void MeshFunctionSnapshot<T>::restore_vector_function(
    const typename mesh::Tria<MeshConfig> &mesh_topology,
    const typename result_of::dof_map_t<MeshConfig> &dof_handler, const SnapshotData &snapshot,
    interpolation::VectorMeshFunctionBase<MFT> &mesh_function)
{
  mesh::StdRegion std_region;
  Uint nb_entries = 0;

  for (Uint i = 0; i < snapshot.cell_paths.size(); ++i)
  {
    const mesh::CellTopologyView<MeshConfig> tcell =
        mesh_topology.path_leaf(snapshot.cell_paths[i]);
    if (tcell.status() == mesh::EntityStatus::Active)
    {
      const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(tcell.active_idx()));
      nb_entries += cell.nb_vert();
      // nb_entries += snapshot.offsets[i + 1] - snapshot.offsets[i];
    }
    else
    {
      const std::vector<mesh::CellTopologyView<MeshConfig>> cell_children = tcell.children();

      for (Uint c = 0; c < cell_children.size(); ++c)
      {
        const mesh::MeshEntity child_cell =
            dof_handler.active_cell(mesh::ActiveIdx(cell_children[c].active_idx()));
        nb_entries += child_cell.nb_vert();
      }
    }
  }

  MFT &function_wrapped = mesh_function.wrapped_type();
  function_wrapped.resize(snapshot.nb_fields, nb_entries);

  ElementAdaptInterpolator local_h_interpolator;
  detail::LocalElemPInterpolation<T> local_p_interpolator;

  math::DenseDMat<T> parent_values;

  for (Uint i = 0; i < snapshot.cell_paths.size(); ++i)
  {
    const mesh::CellTopologyView<MeshConfig> tcell =
        mesh_topology.path_leaf(snapshot.cell_paths[i]);

    // CASE I: if the element remains active, then it was not adapted, or it
    // was p-adapted
    if (tcell.status() == mesh::EntityStatus::Active)
    {
      const mesh::MeshEntity cell = dof_handler.active_cell(mesh::ActiveIdx(tcell.active_idx()));

      const mesh::PointSetTag old_cell_tag = snapshot.cell_type[i].get().pt_set_id();

      // If the tag of the cell is unchanged, then the cell was not
      // adapted (in any way)
      if (old_cell_tag == cell.pt_set_id())
      {
        const Uint nb_vert_in_cell = cell.nb_vert();

        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          typename MFT::entry_type node_value = function_wrapped.value(cell.vertex(v));

          for (Uint f = 0; f < snapshot.nb_fields; ++f)
          {
            node_value[f] = snapshot.values[snapshot.offsets[i] + nb_vert_in_cell * f + v];
          }
        }
      }
      // Else it was p-adapted
      {
        // Put the cached values (before adaptation) into a matrix
        const Uint nb_parent_nodes = snapshot.cell_type[i].get().nb_nodes();
        parent_values.resize(nb_parent_nodes, snapshot.nb_fields);

        const Uint cached_value_idx = snapshot.offsets[i];

        for (Uint n = 0; n < nb_parent_nodes; ++n)
        {
          for (Uint f = 0; f < snapshot.nb_fields; ++f)
          {
            parent_values(n, f) = snapshot.values[cached_value_idx + nb_parent_nodes * f + n];
          }
        }

        local_p_interpolator.interpolate(old_cell_tag, cell.pt_set_id(), parent_values);
        const math::DenseConstMatView<T> interpolated_values =
            local_p_interpolator.interpolated_values();

        for (Uint v = 0; v < cell.nb_vert(); ++v)
        {
          typename MFT::entry_type node_value = function_wrapped.value(cell.vertex(v));

          for (Uint f = 0; f < snapshot.nb_fields; ++f)
          {
            node_value[f] = interpolated_values(v, f);
          }
        }
      }
    }
    // CASE II: if the element is no longer active, then it was h-adapted
    else
    {
      // Put the cached values (before adaptation) into a matrix
      const Uint nb_parent_nodes = snapshot.cell_type[i].get().nb_nodes();
      parent_values.resize(nb_parent_nodes, snapshot.nb_fields);

      const Uint cached_value_idx = snapshot.offsets[i];

      for (Uint n = 0; n < nb_parent_nodes; ++n)
      {
        for (Uint f = 0; f < snapshot.nb_fields; ++f)
        {
          parent_values(n, f) = snapshot.values[cached_value_idx + nb_parent_nodes * f + n];
        }
      }

      // Interpolate the (old) parent values on children
      local_h_interpolator.compute_child_values(snapshot.cell_type[i].get().pt_set_id(),
                                                tcell.cell_adapt_op().get().cell_adapt_op_tag(),
                                                parent_values);

      // Copy the values into the function
      const std::vector<mesh::CellTopologyView<MeshConfig>> cell_children = tcell.children();

      for (Uint c = 0; c < cell_children.size(); ++c)
      {
        const mesh::MeshEntity child_cell =
            dof_handler.active_cell(mesh::ActiveIdx(cell_children[c].active_idx()));
        const math::DenseConstMatView<T> child_values =
            local_h_interpolator.interpolated_child_values(c);
        for (Uint v = 0; v < child_cell.nb_vert(); ++v)
        {
          typename MFT::entry_type node_value = function_wrapped.value(child_cell.vertex(v));

          for (Uint f = 0; f < snapshot.nb_fields; ++f)
          {
            node_value[f] = child_values(v, f);
          }
        }
      } // Loop over cell children
    }   // else
  }     // loop over cell_paths
}

// ----------------------------------------------------------------------------

} // Namespace interpolation

} // Namespace pdekit

#endif
