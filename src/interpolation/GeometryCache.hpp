#ifndef PDEKIT_Interpolation_Geometry_Cache_hpp
#define PDEKIT_Interpolation_Geometry_Cache_hpp

#include <unordered_map>

#include "interpolation/GeometryCacheBase.hpp"
#include "math/DenseConstMatView.hpp"
#include "mesh/adaptation/LocalInterpolator.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode = CacheInsertManual>
class GeometryCache
    : public GeometryCacheBase<GeometryCache<GDIM, CacheInsertMode>, CacheInsertMode>
{
  public:
  enum
  {
    GeoDim = GDIM
  };

  /// Default constructor
  GeometryCache();

  /// Default constructor
  ~GeometryCache();

  /// Empty the data in the buffer
  void flush();

  /// Clear all internal data
  void clear();

  /// Return nb. values that have been pushed to buffer so far
  Uint nb_values_in_buffer() const;

  /// Return the max number of values that the buffer can take
  Uint capacity() const;

  /// Get values corresponding to one element type
  const math::DenseConstMatView<Real> buffer_data(const mesh::DiscreteElemKey key) const;

  /// Get the values corresponding to one cell
  const math::DenseConstMatView<Real> cell_values(const Uint idx) const;

  /// Get the cell type of i-th cell
  const mesh::DiscreteElemKey key(const Uint idx) const;

  /// Position of i-th cell in its data block
  Uint position_in_block(const Uint idx) const;

  /// Print all types that can be pushed to this cache
  void print_types() const;

  /// Print the contents of geometry cache (for debugging)
  void print_contents() const;

  private:
  /// FRIENDS
  friend class GeometryCacheBase<GeometryCache<GDIM, CacheInsertMode>, CacheInsertMode>;

  /// METHODS
  /// Set the size of buffer
  template <typename DiscreteElemKeyIterator>
  void allocate_manual(const DiscreteElemKeyIterator keys_begin,
                       const DiscreteElemKeyIterator keys_end, Uint const nb_blocks);

  /// Insert  block of coordinates to buffer
  void push_back_to_buffer_manual(const mesh::DofCoordinates<GDIM> &cell_coords,
                                  mesh::DiscreteElemKey const key);

  /// Insert  block of coordinates to buffer
  void push_back_to_buffer_manual(const mesh::CellGeometry<GDIM> &cell_coords,
                                  mesh::DiscreteElemKey const key);

  /// Insert  block of coordinates to buffer
  void push_back_to_buffer_manual(const math::DenseConstMatView<Real> &cell_coords,
                                  mesh::DiscreteElemKey const key);

  /// Insert  block of coordinates to buffer
  void push_back_to_buffer_manual_and_interpolate(const mesh::CellGeometry<GDIM> &cell_coords,
                                                  mesh::DiscreteElemKey const coords_key,
                                                  mesh::DiscreteElemKey const interp_key);

  /// Record the order, in which values with different DiscreteElementKeys are
  /// inserted
  /// @note: DiscreteElemKeyInserter must have 3 methods:
  ///        => void init(), which initializes it
  ///        => DiscreteElemKey make_next(), which returns a new
  ///        DiscreteElemKey
  ///        => bool reached_end(), which is true when the inserter finished
  ///        generating
  ///           the insertion parttern
  template <typename DiscreteElemKeyInserter>
  void record_insertion_pattern_automatic(DiscreteElemKeyInserter &inserter, Uint const nb_blocks);

  /// TYPES
  struct BufferData
  {
    /// Type of elements in this buffer data
    mesh::DiscreteElemKey key;

    /// How many elements have been filled
    Uint nb_elem_filled;

    /// Number of dof per elem
    Uint nb_dof_per_elem;

    /// The actual cached values
    math::DenseDMat<Real> m_values;
  };

  using data_map_type =
      std::unordered_map<mesh::DiscreteElemKey, BufferData, mesh::DiscreteElemKeyHash>;

  /// Maximum number of blocks (per element type) that can be stored
  Uint m_max_nb_blocks;

  /// Map of all data arrays corresponding to all element types
  data_map_type m_data_map;

  /// Linear index of metric data, consists of a vector of pairs
  /// [PointSetTag,Uint]
  /// Value m_mdata_index[i] is a pair [FIRST,SECOND] and corresponds to
  /// buffer data which can be found in m_data_map.std_region_data(FIRST),
  /// cell block on position [SECOND] in this particular MetricData
  std::vector<std::pair<typename data_map_type::const_iterator, Uint>> m_data_index;

  /// Global insertion position
  Uint m_next_insert_type_pos;

  /// Interpolator for transferring the coordinates from one point set to another
  /// (provided element shape remains the same)
  mesh::adapt::LocalInterpolator m_loc_interp;
};

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
GeometryCache<GDIM, CacheInsertMode>::GeometryCache()
    : m_max_nb_blocks(0U), m_next_insert_type_pos(0)
{
  m_data_map.clear();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
GeometryCache<GDIM, CacheInsertMode>::~GeometryCache()
{
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::flush()
{
  for (typename data_map_type::iterator it = m_data_map.begin(); it != m_data_map.end(); ++it)
  {
    (*it).second.nb_elem_filled = 0;
  }
  m_data_index.resize(0);
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::clear()
{
  m_max_nb_blocks = 0;
  m_data_map.clear();
  m_data_index.clear();
  m_next_insert_type_pos = 0;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
Uint GeometryCache<GDIM, CacheInsertMode>::nb_values_in_buffer() const
{
  return m_data_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
Uint GeometryCache<GDIM, CacheInsertMode>::capacity() const
{
  return m_max_nb_blocks;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
inline const math::DenseConstMatView<Real> GeometryCache<GDIM, CacheInsertMode>::buffer_data(
    const mesh::DiscreteElemKey key) const
{
  typename data_map_type::const_iterator buffer_data_it = m_data_map.find(key);
  math::DenseDMat<Real> const &mat                      = (*buffer_data_it).second.m_values;

  math::DenseConstMatView<Real> block(&mat(0, 0), mat.cols(), mat.rows(),
                                      (*buffer_data_it).second.nb_elem_filled * GDIM);
  return block;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
inline const math::DenseConstMatView<Real> GeometryCache<GDIM, CacheInsertMode>::cell_values(
    const Uint idx) const
{
  typename data_map_type::const_iterator buffer_data_it = m_data_index[idx].first;
  math::DenseDMat<Real> const &mat                      = (*buffer_data_it).second.m_values;

  math::DenseConstMatView<Real> cell_values(&mat(0, m_data_index[idx].second * GDIM), mat.cols(),
                                            mat.rows(), GDIM);

  // math::MatrixBlock<Real> cell_values(
  //    mat.block(0, 0, mat.rows(), GDIM));

  return cell_values;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
inline const mesh::DiscreteElemKey GeometryCache<GDIM, CacheInsertMode>::key(const Uint idx) const
{
  const typename data_map_type::const_iterator buffer_data = m_data_index[idx].first;
  return (*buffer_data).first;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
inline Uint GeometryCache<GDIM, CacheInsertMode>::position_in_block(const Uint idx) const
{
  return m_data_index[idx].second;
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::print_types() const
{
  for (typename data_map_type::const_iterator it = m_data_map.cbegin(); it != m_data_map.cend();
       ++it)
  {
    std::cout << "[" << (*it).first << "], fill status " << (*it).second.nb_elem_filled << "/"
              << m_max_nb_blocks << " items" << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::print_contents() const
{
  for (typename data_map_type::const_iterator it = m_data_map.cbegin(); it != m_data_map.cend();
       ++it)
  {
    std::cout << "[" << (*it).first << "], fill status " << (*it).second.nb_elem_filled << "/"
              << m_max_nb_blocks << " items" << std::endl;
    std::cout << (*it).second.m_values << std::endl << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
template <typename DiscreteElemKeyIterator>
void GeometryCache<GDIM, CacheInsertMode>::allocate_manual(const DiscreteElemKeyIterator keys_begin,
                                                           const DiscreteElemKeyIterator keys_end,
                                                           Uint const nb_blocks)

{
  m_max_nb_blocks = nb_blocks;
  m_data_map.clear();
  m_next_insert_type_pos = 0;

  DiscreteElemKeyIterator key_iter = keys_begin;
  Uint n_elem_types                = 0;

  for (; key_iter != keys_end; ++key_iter)
  {
    n_elem_types++;
  }

  m_data_index.reserve(n_elem_types * nb_blocks);
  m_data_index.resize(0);

  mesh::StdRegion ref_element;

  key_iter = keys_begin;

  for (; key_iter != keys_end; ++key_iter)
  {
    const mesh::PointSetTag std_region_tag = (*key_iter).support().std_region_tag();
    // const Uint local_id = elem_iter.key_value().local_id;

    ref_element.change_type(std_region_tag);

    // common::PtrHandle<BufferData> buffer_data =
    //     m_data_map.create(mesh::DefaultStdRegMapKey(std_region_tag,
    //     local_id, P0));

    BufferData buffer_data;

    const Uint nb_dof_per_elem = ref_element.get().nb_nodes();

    buffer_data.key             = *key_iter;
    buffer_data.nb_dof_per_elem = nb_dof_per_elem;
    buffer_data.nb_elem_filled  = 0U;
    buffer_data.m_values.resize(nb_dof_per_elem, GDIM * m_max_nb_blocks);

    m_data_map.insert(std::move(std::make_pair(*key_iter, buffer_data)));
  }
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::push_back_to_buffer_manual(
    const mesh::DofCoordinates<GDIM> &cell_coords, mesh::DiscreteElemKey const key)
{
  using point_coord_type = typename mesh::DofCoordinates<GDIM>::node_coords;

  typename data_map_type::iterator buffer_data_it  = m_data_map.find(key);
  typename data_map_type::mapped_type &buffer_data = (*buffer_data_it).second;
  const Uint offset                                = buffer_data.nb_elem_filled * GDIM;

  math::DenseDMat<Real> &active_buffer = buffer_data.m_values;

  /// @TODO: this loop should not go till cell_coords.size(), but only to the
  /// number of nodes corresponding to the number of nodes of standard region
  /// of type 'cell_type' Consider what would happen if cell_coords correspond
  /// to a volume element, but cell_type corresponds to a standard region type
  /// of one of its faces!
  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    point_coord_type const node = cell_coords.c(n);

    for (Uint d = 0; d < GDIM; ++d)
    {
      active_buffer(n, offset + d) = node[d];
    }
  }

  m_data_index.push_back(std::pair<typename data_map_type::const_iterator, Uint>(
      buffer_data_it, buffer_data.nb_elem_filled));

  buffer_data.nb_elem_filled++;
  m_next_insert_type_pos = (m_next_insert_type_pos + 1) % m_data_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::push_back_to_buffer_manual(
    const mesh::CellGeometry<GDIM> &cell_coords, mesh::DiscreteElemKey const key)
{
  using point_coord_type = typename mesh::CellGeometry<GDIM>::node_view_t;

  typename data_map_type::iterator buffer_data_it  = m_data_map.find(key);
  typename data_map_type::mapped_type &buffer_data = (*buffer_data_it).second;

  const Uint offset = buffer_data.nb_elem_filled * GDIM;

  math::DenseDMat<Real> &active_buffer = buffer_data.m_values;

  /// @TODO: this loop should not go till cell_coords.size(), but only to the
  /// number of nodes corresponding to the number of nodes of standard region
  /// of type 'cell_type' Consider what would happen if cell_coords correspond
  /// to a volume element, but cell_type corresponds to a standard region type
  /// of one of its faces!
  for (Uint n = 0; n < cell_coords.size(); ++n)
  {
    point_coord_type const node = cell_coords.const_node_view(n);

    for (Uint d = 0; d < GDIM; ++d)
    {
      active_buffer(n, offset + d) = node[d];
    }
  }

  m_data_index.push_back(std::pair<typename data_map_type::const_iterator, Uint>(
      buffer_data_it, buffer_data.nb_elem_filled));

  buffer_data.nb_elem_filled++;
  m_next_insert_type_pos = (m_next_insert_type_pos + 1) % m_data_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::push_back_to_buffer_manual(
    const math::DenseConstMatView<Real> &cell_coords, mesh::DiscreteElemKey const key)
{
  typename data_map_type::iterator buffer_data_it  = m_data_map.find(key);
  typename data_map_type::mapped_type &buffer_data = (*buffer_data_it).second;
  const Uint offset                                = buffer_data.nb_elem_filled * GDIM;

  math::DenseDMat<Real> &active_buffer = buffer_data.m_values;

  /// @TODO: this loop should not go till cell_coords.size(), but only to the
  /// number of nodes corresponding to the number of nodes of standard region
  /// of type 'cell_type' Consider what would happen if cell_coords correspond
  /// to a volume element, but cell_type corresponds to a standard region type
  /// of one of its faces!
  for (Uint n = 0; n < cell_coords.rows(); ++n)
  {
    for (Uint d = 0; d < GDIM; ++d)
    {
      active_buffer(n, offset + d) = cell_coords(n, d);
    }
  }

  m_data_index.push_back(std::pair<typename data_map_type::const_iterator, Uint>(
      buffer_data_it, buffer_data.nb_elem_filled));

  buffer_data.nb_elem_filled++;
  m_next_insert_type_pos = (m_next_insert_type_pos + 1) % m_data_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
void GeometryCache<GDIM, CacheInsertMode>::push_back_to_buffer_manual_and_interpolate(
    const mesh::CellGeometry<GDIM> &cell_coords, mesh::DiscreteElemKey const coords_key,
    mesh::DiscreteElemKey const interp_key)
{
  typename data_map_type::iterator buffer_data_it  = m_data_map.find(interp_key);
  typename data_map_type::mapped_type &buffer_data = (*buffer_data_it).second;

  const Uint offset = buffer_data.nb_elem_filled * GDIM;

  math::DenseDMat<Real> &active_buffer = buffer_data.m_values;

  const math::DenseConstMatView<Real> interp_coords =
      m_loc_interp.transfer_coords(coords_key, interp_key, cell_coords);

  for (Uint n = 0; n < interp_coords.rows(); ++n)
  {
    for (Uint d = 0; d < GDIM; ++d)
    {
      active_buffer(n, offset + d) = interp_coords(n, d);
    }
  }

  m_data_index.push_back(std::pair<typename data_map_type::const_iterator, Uint>(
      buffer_data_it, buffer_data.nb_elem_filled));

  buffer_data.nb_elem_filled++;
  m_next_insert_type_pos = (m_next_insert_type_pos + 1) % m_data_index.size();
}

// ----------------------------------------------------------------------------

template <Uint GDIM, typename CacheInsertMode>
template <typename DiscreteElemKeyInserter>
void GeometryCache<GDIM, CacheInsertMode>::record_insertion_pattern_automatic(
    DiscreteElemKeyInserter &inserter, Uint const nb_blocks)
{
  m_max_nb_blocks = nb_blocks;
  m_data_map.clear();
  m_next_insert_type_pos = 0;

  mesh::StdRegion ref_element;

  Uint n_inserted_values = 0;

  // First run: allocated data
  inserter.init();
  while (!inserter.reached_end())
  {
    const mesh::DiscreteElemKey key     = inserter.make_next();
    typename data_map_type::iterator it = m_data_map.find(key);

    if (it == m_data_map.end())
    {
      const mesh::PointSetTag std_region_tag = key.support().std_region_tag();
      ref_element.change_type(std_region_tag);
      const Uint nb_dof_per_elem = ref_element.get().nb_nodes();

      BufferData buffer_data;
      buffer_data.key             = key;
      buffer_data.nb_dof_per_elem = nb_dof_per_elem;
      buffer_data.nb_elem_filled  = 0;
      buffer_data.m_values.resize(nb_dof_per_elem, GDIM * m_max_nb_blocks);
      m_data_map.insert(std::move(std::make_pair(key, buffer_data)));
    }
    n_inserted_values++;
  }

  const Uint n_elem_types = m_data_map.size();

  m_data_index.reserve(n_elem_types * nb_blocks);
  m_data_index.resize(0);

  // Second run: record order in which values are inserted
  inserter.init();
  while (!inserter.reached_end())
  {
    const mesh::DiscreteElemKey key                 = inserter.make_next();
    typename data_map_type::iterator buffer_data_it = m_data_map.find(key);

    typename data_map_type::mapped_type &buffer_data = (*buffer_data_it).second;

    m_data_index.push_back(std::pair<typename data_map_type::const_iterator, Uint>(
        buffer_data_it, buffer_data.nb_elem_filled));

    buffer_data.nb_elem_filled++;
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit

#endif
