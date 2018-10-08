#ifndef PDEKIT_Interpolation_Solution_Cache_hpp
#define PDEKIT_Interpolation_Solution_Cache_hpp

#include "common/DataMap.hpp"
#include "common/Meta.hpp"
#include "common/Range1D.hpp"
#include "interpolation/FEValues.hpp"
#include "math/DenseConstMatView.hpp"
#include "math/DenseConstVecView.hpp"
#include "mesh/EntityDofRealign.hpp"
#include "mesh/MeshEntity.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace interpolation
{

template <typename T>
class ScalarMeshFunction;
template <typename T>
class VectorMeshFunction;

class SolutionCache
{
  public:
  /// Default constructor
  SolutionCache();

  /// Default constructor
  ~SolutionCache();

  /// Set the size of buffer
  /// @param elem_types ... map of FiniteElement values associated with each
  /// element type
  /// @param nb_blocks ... maximum number of elements (of each type) that the
  /// cache can hold
  /// @param nb_fields ... number of components of the solution field when
  /// solving a coupled
  ///                      (i.e. not scalar) problem
  template <typename FeValsIterator>
  void allocate(const FeValsIterator fe_vals_begin, const FeValsIterator fe_vals_end,
                Uint const nb_blocks, const Uint nb_fields);

  /// Copy the cache to another one. This does not allocate data (has to be
  /// done previously) The only thing copied is the contents of the containers
  /// @param cache_in ... solution cache to be copied
  /// @param cache_out ... target cache of the copying process
  static void copy_contents(const SolutionCache &cache_in, SolutionCache &cache_out);

  /// Insert generic data to buffer
  /// @param cell ... cell to be inserted - has information about element
  /// type, node numbers etc.
  ///                 The node indexes determine where in the next input
  ///                 parameter 'data' are the values to be copied
  /// @param data ... container from which the data is being inserted
  /// @param cell_type ... element type associated with this data
  template <typename T>
  void push_back_to_buffer(mesh::MeshEntity const &cell, ScalarMeshFunction<T> const &data,
                           mesh::PointSetTagExt const cell_type);

  /// Insert generic data to buffer
  /// @param cell ... cell to be inserted - has information about element
  /// type, node numbers etc.
  ///                 The node indexes determine where in the next input
  ///                 parameter 'data' are the values to be copied
  /// @param p ... permutation of the nodes of cell (previous parameter) which
  /// determines the order
  ///              in which the data values are inserted in the cache (there
  ///              is always one solution value/vector associated with one
  ///              node of cell, but the order of cell nodes can be permuted)
  /// @param data ... container from which the data is being inserted
  /// @param cell_type ... element type associated with this data
  template <typename T>
  void push_back_to_buffer(mesh::MeshEntity const &cell, mesh::EntityDofRealign p,
                           ScalarMeshFunction<T> const &data, mesh::PointSetTagExt const cell_type);

  /// Insert generic data to buffer
  /// @param cell ... cell to be inserted - has information about element
  /// type, node numbers etc.
  ///                 The node indexes determine where in the next input
  ///                 parameter 'data' are the values to be copied
  /// @param data ... container from which the data is being inserted
  /// @param cell_type ... element type associated with this data
  template <typename T>
  void push_back_to_buffer(mesh::MeshEntity const &cell, VectorMeshFunction<T> const &data,
                           mesh::PointSetTagExt const cell_type);

  /// Insert generic data to buffer
  /// @param cell ... cell to be inserted - has information about element
  /// type, node numbers etc.
  ///                 The node indexes determine where in the next input
  ///                 parameter 'data' are the values to be copied
  /// @param p ... permutation of the nodes of cell (previous parameter) which
  /// determines the order
  ///              in which the data values are inserted in the cache (there
  ///              is always one solution value/vector associated with one
  ///              node of cell, but the order of cell nodes can be permuted)
  /// @param data ... container from which the data is being inserted
  /// @param cell_type ... element type associated with this data
  template <typename T>
  void push_back_to_buffer(mesh::MeshEntity const &cell, mesh::EntityDofRealign p,
                           VectorMeshFunction<T> const &data, mesh::PointSetTagExt const cell_type);

  template <typename VT, bool TF>
  void push_vec_to_buffer(mesh::MeshEntity const &cell, math::DenseVector<VT, TF> const &data,
                          mesh::PointSetTagExt const cell_type);

  /// Perform perturbation of solution values in one node of each element of
  /// given type (used for implicit methods)
  /// @param std_region ... type of element whose solution values should be
  /// perturbed
  ///                       Solution values associated with other element
  ///                       types than std_region which are currently held in
  ///                       cache will remain unmodified
  /// @param local_node ... node index within ONE element of type std_region
  /// whose associated
  ///                       solution value should be perturbed. Note that the
  ///                       solution value associated with one node can be
  ///                       VECTOR - we have to select its component
  /// @param component ... component of vector associated with local_node in
  /// one element which should
  ///                      be perturbed
  /// @param factor ... perturbation factor
  /*
  void perturb_values(const mesh::PointSetTagExt std_region, const Uint
  local_node, const Uint component, const Real factor = 1.e-7);
  */

  /// Perform perturbation of solution values in one node of each element of
  /// given type (used for implicit methods)
  /// @param elem_range ... range of elements whose solution values should be
  /// perturbed
  ///                       Solution values associated with other element
  ///                       types than std_region which are currently held in
  ///                       cache will remain unmodified
  /// @param local_node ... node index within ONE element of type std_region
  /// whose associated
  ///                       solution value should be perturbed. Note that the
  ///                       solution value associated with one node can be
  ///                       VECTOR - we have to select its component
  /// @param component ... component of vector associated with local_node in
  /// one element which should
  ///                      be perturbed
  /// @param factor ... perturbation factor
  void perturb_values(const common::Range1D<Uint> elem_range, const Uint local_node,
                      const Uint component, const Real factor = 1.e-7);

  /// Reset the cache values so that nothing is perturbed
  /// @param std_region ... type of element whose data should be reset to 'no
  /// perturbation'
  ///                       Solution data associated with other element types
  ///                       stored in the cache will not be affected by this
  ///                       operation
  // void remove_perturbation(const mesh::PointSetTagExt std_region);

  /// Reset the cache values so that nothing is perturbed
  /// @param elem_range ... range of elements that should be reset to 'no
  /// perturbation'
  ///                       Solution data associated with other elements
  ///                       stored in the cache will not be affected by this
  ///                       operation
  void remove_perturbation(const common::Range1D<Uint> elem_range);

  /// Get the original solution values for given cell type - one value per
  /// cell The returned vector contains values for an array of cells of the
  /// same type
  /*
  const math::DenseConstVecView<Real>
  unperturbed_values(const mesh::PointSetTagExt std_region) const;
  */

  /// Get the original solution values for given cell type - one value per
  /// cell The returned vector contains values for an array of cells of the
  /// same type
  const math::DenseConstVecView<Real> unperturbed_values(
      const common::Range1D<Uint> elem_range) const;

  /// Empty the data in the buffer - set the data so that
  /// no solution values are stored in the cache
  void flush();

  /// Clear all internal data - pre-allocated buffers will be destroyed and
  /// all element types that can be recognized by the cache will be
  /// 'forgotten'
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Maximum number of values that the cache can hold
  Uint capacity() const;

  /// Get values corresponding to one element type
  const math::DenseConstMatView<Real> buffer_data(const mesh::PointSetTagExt std_region) const;

  /// Get the values corresponding to one cell
  const math::DenseConstMatView<Real> cell_values(const Uint idx) const;

  /// Get the cell type of i-th cell
  const mesh::PointSetTagExt std_region_type(const Uint idx) const;

  /// Position of i-th cell in its data block
  Uint position_in_block(const Uint idx) const;

  /// Print all types that can be pushed to this cache
  void print_types() const;

  private:
  /// TYPES

  struct BufferData
  {
    /// Type of elements in this buffer data
    mesh::PointSetTagExt std_region_type;

    /// How many elements have been filled
    Uint nb_elem_filled;

    /// Number of dof per elem
    Uint nb_dof_per_elem;

    /// The actual cached values
    math::DenseDMat<Real> m_values;

    /// Offset indicating where to find the perturbation data
    /// in arrays:
    /// m_perturb_cache
    /// local_perturb_node
    /// local_perturb_component
    Uint m_perturb_data_first;
    Uint m_perturb_data_last;
  };

  /// Maximum number of blocks (per element type) that can be stored
  Uint m_max_nb_blocks;

  /// Number of fields per block
  Uint m_nb_fields;

  /// Map of all data arrays corresponding to all element types
  common::DataMap<mesh::PointSetTagExt, BufferData> m_data_map;

  /// Perturbation - one perturbed value per cell, i.e.
  /// if we store say four-component solution on P2 triangles,
  /// we have blocks of size 6 x 4 values (6 nodes, 4 values per node),
  /// but only ONE component in one selected node can be perturbed at a time
  math::DenseDVec<Real> m_perturb_cache;

  /// Local index of dof in ONE element which is perturbed
  math::DenseDVec<Uint> m_local_perturb_node;

  /// Component of solution in local_perturb_node which is perturbed
  math::DenseDVec<Uint> m_local_perturb_component;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// buffer data which can be found in m_data_map.std_region_data(FIRST),
  /// cell block on position [SECOND] in this particular MetricData
  std::vector<std::pair<common::PtrHandle<BufferData>, Uint>> m_data_index;
};

// ----------------------------------------------------------------------------

template <typename FeValsIterator>
void SolutionCache::allocate(const FeValsIterator fe_vals_begin, const FeValsIterator fe_vals_end,
                             Uint const nb_blocks, const Uint nb_fields)
{
  m_max_nb_blocks = nb_blocks;
  m_nb_fields     = nb_fields;

  FeValsIterator it = fe_vals_begin;
  Uint n_elem_types = 0;

  for (; it != fe_vals_end; ++it)
  {
    n_elem_types++;
  }

  m_data_index.reserve(n_elem_types * nb_blocks);
  m_data_index.resize(0);

  mesh::StdRegion ref_element;

  Uint perturb_data_offset = 0;

  it = fe_vals_begin;

  for (; it != fe_vals_end; ++it)
  {
    const mesh::PointSetTag std_region_tag = it.key_value().std_region_tag();
    // const Uint local_id = it.key_value().local_id;

    ref_element.change_type(std_region_tag);

    common::PtrHandle<BufferData> buffer_data = m_data_map.create(it.key_value());

    const Uint nb_dof_per_elem = ref_element.get().nb_nodes();

    (*buffer_data).std_region_type = it.key_value();
    (*buffer_data).nb_elem_filled  = 0U;
    (*buffer_data).nb_dof_per_elem = nb_dof_per_elem;
    (*buffer_data).m_values.resize(nb_dof_per_elem, m_nb_fields * m_max_nb_blocks);
    (*buffer_data).m_perturb_data_first = perturb_data_offset;
    (*buffer_data).m_perturb_data_last  = perturb_data_offset + m_max_nb_blocks - 1;
    perturb_data_offset += m_max_nb_blocks;
  }
  m_perturb_cache.resize(m_max_nb_blocks * n_elem_types);
  m_local_perturb_node.resize(m_max_nb_blocks * n_elem_types);
  m_local_perturb_component.resize(m_max_nb_blocks * n_elem_types);
}

// ----------------------------------------------------------------------------

template <typename T>
void SolutionCache::push_back_to_buffer(mesh::MeshEntity const &cell,
                                        ScalarMeshFunction<T> const &data,
                                        mesh::PointSetTagExt const cell_type)
{
  common::PtrHandle<BufferData> buffer_data = m_data_map.std_region_data(cell_type);
  const Uint offset                         = (*buffer_data).nb_elem_filled * m_nb_fields;

  math::DenseDMat<Real> &active_buffer = (*buffer_data).m_values;

  for (Uint n = 0; n < cell.nb_vert(); ++n)
  {
    active_buffer(n, offset + 0) = data[cell.vertex(n)];
  }

  m_data_index.push_back(
      std::pair<common::PtrHandle<BufferData>, Uint>(buffer_data, (*buffer_data).nb_elem_filled));
  (*buffer_data).nb_elem_filled++;
}

// ----------------------------------------------------------------------------

template <typename T>
void SolutionCache::push_back_to_buffer(mesh::MeshEntity const &cell, mesh::EntityDofRealign p,
                                        ScalarMeshFunction<T> const &data,
                                        mesh::PointSetTagExt const cell_type)
{
  common::PtrHandle<BufferData> buffer_data = m_data_map.std_region_data(cell_type);
  const Uint offset                         = (*buffer_data).nb_elem_filled * m_nb_fields;

  math::DenseDMat<Real> &active_buffer = (*buffer_data).m_values;

  for (Uint n = 0; n < cell.nb_vert(); ++n)
  {
    active_buffer(n, offset + 0) = data[cell.vertex(p.get().vertex(n))];
  }

  m_data_index.push_back(
      std::pair<common::PtrHandle<BufferData>, Uint>(buffer_data, (*buffer_data).nb_elem_filled));
  (*buffer_data).nb_elem_filled++;
}

// ----------------------------------------------------------------------------

template <typename T>
void SolutionCache::push_back_to_buffer(mesh::MeshEntity const &cell,
                                        VectorMeshFunction<T> const &data,
                                        mesh::PointSetTagExt const cell_type)
{
  // std::cout << "Filling buffer with cell " << cell << std::endl;
  common::PtrHandle<BufferData> buffer_data = m_data_map.std_region_data(cell_type);
  const Uint offset                         = (*buffer_data).nb_elem_filled * m_nb_fields;

  math::DenseDMat<Real> &active_buffer = (*buffer_data).m_values;

  for (Uint n = 0; n < cell.nb_vert(); ++n)
  {
    typename VectorMeshFunction<T>::const_entry_type nodal_value = data.const_value(cell.vertex(n));

    for (Uint f = 0; f < m_nb_fields; ++f)
    {
      active_buffer(n, offset + f) = nodal_value[f];
    }
  }

  m_data_index.push_back(
      std::pair<common::PtrHandle<BufferData>, Uint>(buffer_data, (*buffer_data).nb_elem_filled));
  (*buffer_data).nb_elem_filled++;
}

// ----------------------------------------------------------------------------

template <typename T>
void SolutionCache::push_back_to_buffer(mesh::MeshEntity const &cell, mesh::EntityDofRealign p,
                                        VectorMeshFunction<T> const &data,
                                        mesh::PointSetTagExt const cell_type)
{
  // std::cout << "Filling buffer with cell " << cell << std::endl;

  common::PtrHandle<BufferData> buffer_data = m_data_map.std_region_data(cell_type);
  const Uint offset                         = (*buffer_data).nb_elem_filled * m_nb_fields;

  math::DenseDMat<Real> &active_buffer = (*buffer_data).m_values;

  for (Uint n = 0; n < cell.nb_vert(); ++n)
  {
    typename VectorMeshFunction<T>::const_entry_type data_value =
        data.const_value(cell.vertex(p.get().vertex(n)));

    for (Uint f = 0; f < m_nb_fields; ++f)
    {
      active_buffer(n, offset + f) = data_value[f];
    }
  }

  m_data_index.push_back(
      std::pair<common::PtrHandle<BufferData>, Uint>(buffer_data, (*buffer_data).nb_elem_filled));
  (*buffer_data).nb_elem_filled++;
}

// ----------------------------------------------------------------------------

template <typename VT, bool TF>
void SolutionCache::push_vec_to_buffer(mesh::MeshEntity const &cell,
                                       math::DenseVector<VT, TF> const &data,
                                       mesh::PointSetTagExt const cell_type)
{
  common::PtrHandle<BufferData> buffer_data = m_data_map.std_region_data(cell_type);
  const Uint offset                         = (*buffer_data).nb_elem_filled * m_nb_fields;

  const VT &vec_data = data.wrapped_type();

  math::DenseDMat<Real> &active_buffer = (*buffer_data).m_values;

  for (Uint n = 0; n < cell.nb_vert(); ++n)
  {
    for (Uint f = 0; f < m_nb_fields; ++f)
    {
      active_buffer(n, offset + f) = vec_data[n * m_nb_fields + f];
    }
  }

  m_data_index.push_back(
      std::pair<common::PtrHandle<BufferData>, Uint>(buffer_data, (*buffer_data).nb_elem_filled));
  (*buffer_data).nb_elem_filled++;
}

// ----------------------------------------------------------------------------

const math::DenseConstMatView<Real> inline SolutionCache::buffer_data(
    const mesh::PointSetTagExt std_region) const
{
  common::PtrHandle<BufferData const> buffer_data = m_data_map.std_region_data(std_region);
  math::DenseDMat<Real> const &mat                = (*buffer_data).m_values;

  math::DenseConstMatView<Real> block(&mat(0, 0), mat.cols(), mat.rows(),
                                      (*buffer_data).nb_elem_filled * m_nb_fields);
  return block;
}

// ----------------------------------------------------------------------------

inline const math::DenseConstMatView<Real> SolutionCache::cell_values(const Uint idx) const
{
  common::PtrHandle<BufferData const> buffer_data = m_data_index[idx].first;
  math::DenseDMat<Real> const &mat                = (*buffer_data).m_values;

  math::DenseConstMatView<Real> cell_values(&mat(0, m_data_index[idx].second * m_nb_fields),
                                            mat.cols(), mat.rows(), m_nb_fields);

  /*
  math::MatrixBlock<Real>
  cell_values(mat.block(0, m_data_index[idx].second * m_nb_fields,
  mat.rows(), m_nb_fields));
  */

  return cell_values;
}

// ----------------------------------------------------------------------------

inline const mesh::PointSetTagExt SolutionCache::std_region_type(const Uint idx) const
{
  common::PtrHandle<BufferData const> buffer_data = m_data_index[idx].first;
  return (*buffer_data).std_region_type;
}

// ----------------------------------------------------------------------------

inline Uint SolutionCache::position_in_block(const Uint idx) const
{
  return m_data_index[idx].second;
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
