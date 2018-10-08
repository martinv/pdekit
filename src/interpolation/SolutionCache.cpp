#include "interpolation/SolutionCache.hpp"
#include "math/DenseVecView.hpp"
#include "mesh/std_region/StdRegion.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

SolutionCache::SolutionCache() : m_max_nb_blocks(0U), m_nb_fields(0U)
{
  m_data_map.clear();
}

// ----------------------------------------------------------------------------

SolutionCache::~SolutionCache()
{
}

// ----------------------------------------------------------------------------

void SolutionCache::copy_contents(const SolutionCache &cache_in, SolutionCache &cache_out)
{
  using data_map_type = common::DataMap<mesh::PointSetTagExt, BufferData>;

  for (data_map_type::const_iterator buff_it_in = cache_in.m_data_map.cbegin();
       buff_it_in != cache_in.m_data_map.cend(); ++buff_it_in)
  {
    BufferData const &data_in = (*buff_it_in.data_ptr());
    BufferData &data_out      = (*cache_out.m_data_map.std_region_data(buff_it_in.key_value()));

    data_out.std_region_type = data_in.std_region_type;
    data_out.nb_elem_filled  = data_in.nb_elem_filled;
    data_out.nb_dof_per_elem = data_in.nb_dof_per_elem;

    if ((data_out.m_values.rows() != data_in.m_values.rows()) ||
        (data_out.m_values.cols() != data_in.m_values.cols()))
    {
      data_out.m_values.resize(data_in.m_values.rows(), data_in.m_values.cols());
    }
    // Copy the actual data matrices
    data_out.m_values = data_in.m_values;

    data_out.m_perturb_data_first = data_in.m_perturb_data_first;
    data_out.m_perturb_data_last  = data_in.m_perturb_data_last;
  } // Loop over available element types

  if (cache_out.m_perturb_cache.size() != cache_in.m_perturb_cache.size())
  {
    cache_out.m_perturb_cache.resize(cache_in.m_perturb_cache.size());
  }

  if (cache_out.m_local_perturb_node.size() != cache_in.m_local_perturb_node.size())
  {
    cache_out.m_local_perturb_node.resize(cache_in.m_local_perturb_node.size());
  }

  if (cache_out.m_local_perturb_component.size() != cache_in.m_local_perturb_component.size())
  {
    cache_out.m_local_perturb_component.resize(cache_in.m_local_perturb_component.size());
  }

  cache_out.m_perturb_cache.resize(cache_in.m_perturb_cache.size());
  cache_out.m_local_perturb_node.resize(cache_in.m_local_perturb_node.size());
  cache_out.m_local_perturb_component.resize(cache_in.m_local_perturb_component.size());

  cache_out.m_perturb_cache           = cache_in.m_perturb_cache;
  cache_out.m_local_perturb_node      = cache_in.m_local_perturb_node;
  cache_out.m_local_perturb_component = cache_in.m_local_perturb_component;

  cache_out.m_data_index.resize(cache_in.m_data_index.size());
  for (Uint i = 0; i < cache_in.m_data_index.size(); ++i)
  {
    const mesh::PointSetTagExt key_in = (*cache_in.m_data_index[i].first).std_region_type;
    const Uint index_in               = cache_in.m_data_index[i].second;

    common::PtrHandle<BufferData> const buffer_ptr = cache_out.m_data_map.std_region_data(key_in);
    cache_out.m_data_index.push_back(
        std::pair<common::PtrHandle<BufferData>, Uint>(buffer_ptr, index_in));
  }
}

// ----------------------------------------------------------------------------

/*
void SolutionCache::perturb_values(const mesh::PointSetTagExt std_region, const
Uint local_node, const Uint component, const Real factor)
{
  BufferData &buffer_data = (*m_data_map.std_region_data(std_region));
  const Uint perturb_data_offset = buffer_data.m_perturb_data_first;

  // Shape of buffer_data.m_values: [ Nb. of shape functions x (Nb. elem * Nb.
of fields)] for (Uint e = 0; e < buffer_data.nb_elem_filled; ++e)
  {
    const Real u = buffer_data.m_values(local_node, e * m_nb_fields +
component); const Real sign_u = u < 0.0 ? -1.0 : 1.0; const Real u_perturb =
factor * sign_u * std::max(std::abs(u), 0.01);

    // Move the original value to buffer_data.m_perturb_cache[e];
    buffer_data.m_values(local_node, e * m_nb_fields + component) += u_perturb;

    m_perturb_cache[perturb_data_offset + e] = u;
    m_local_perturb_node[perturb_data_offset + e] = local_node;
    m_local_perturb_component[perturb_data_offset + e] = component;
  }
}
*/

// ----------------------------------------------------------------------------

void SolutionCache::perturb_values(const common::Range1D<Uint> elem_range, const Uint local_node,
                                   const Uint component, const Real factor)
{
  for (Uint eidx = elem_range.lbound(); eidx <= elem_range.ubound(); ++eidx)
  {
    BufferData &buffer_data     = *(m_data_index[eidx].first);
    const Uint idx_in_buff_data = m_data_index[eidx].second;

    const Real u = buffer_data.m_values(local_node, idx_in_buff_data * m_nb_fields + component);
    const Real sign_u    = u < 0.0 ? -1.0 : 1.0;
    const Real u_perturb = factor * sign_u * std::max(std::abs(u), 0.01);

    // Move the original value to buffer_data.m_perturb_cache[e];
    buffer_data.m_values(local_node, idx_in_buff_data * m_nb_fields + component) += u_perturb;

    m_perturb_cache[eidx]           = u;
    m_local_perturb_node[eidx]      = local_node;
    m_local_perturb_component[eidx] = component;
  }
}

// ----------------------------------------------------------------------------

/*
void SolutionCache::remove_perturbation(const mesh::PointSetTagExt std_region)
{
  BufferData &buffer_data = (*m_data_map.std_region_data(std_region));
  // const Uint perturb_data_offset = buffer_data.m_perturb_data_first;

  math::DenseVecView<Real>
data_perturb_cache(&m_perturb_cache[buffer_data.m_perturb_data_first],
                                              &m_perturb_cache[buffer_data.m_perturb_data_last]);

  math::DenseVecView<Uint> data_local_perturb_node(
      &m_local_perturb_node[buffer_data.m_perturb_data_first],
      &m_local_perturb_node[buffer_data.m_perturb_data_last]);

  math::DenseVecView<Uint> data_local_perturb_component(
      &m_local_perturb_component[buffer_data.m_perturb_data_first],
      &m_local_perturb_component[buffer_data.m_perturb_data_last]);

  // Shape of buffer_data.m_values: [ Nb. of shape functions x (Nb. elem * Nb.
of fields)] for (Uint e = 0; e < buffer_data.nb_elem_filled; ++e)
  {
    buffer_data.m_values(data_local_perturb_node[e],
                         e * m_nb_fields + data_local_perturb_component[e]) =
data_perturb_cache[e]; data_perturb_cache[e] = 0.0; data_local_perturb_node[e] =
0; data_local_perturb_component[e] = 0;
  }
}
*/

// ----------------------------------------------------------------------------

void SolutionCache::remove_perturbation(const common::Range1D<Uint> elem_range)
{
  for (Uint eidx = elem_range.lbound(); eidx <= elem_range.ubound(); ++eidx)
  {
    BufferData &buffer_data     = *(m_data_index[eidx].first);
    const Uint idx_in_buff_data = m_data_index[eidx].second;

    buffer_data.m_values(m_local_perturb_node[eidx],
                         idx_in_buff_data * m_nb_fields + m_local_perturb_component[eidx]) =
        m_perturb_cache[eidx];
    m_perturb_cache[eidx]           = 0.0;
    m_local_perturb_node[eidx]      = 0;
    m_local_perturb_component[eidx] = 0;
  }
}

// ----------------------------------------------------------------------------

#if 0
const math::DenseConstVecView<Real>
SolutionCache::unperturbed_values(const mesh::PointSetTagExt std_region) const
{
  BufferData const &buffer_data = (*m_data_map.std_region_data(std_region));

  /*
  const math::DenseConstVecView<Real> perturbations(&buffer_data.m_perturb_cache[0],
                                                    buffer_data.m_perturb_cache.size());
  */

  const math::DenseConstVecView<Real> cached_values(
      &m_perturb_cache[buffer_data.m_perturb_data_first],
      &m_perturb_cache[buffer_data.m_perturb_data_last]);
  return cached_values;
}
#endif

// ----------------------------------------------------------------------------

const math::DenseConstVecView<Real> SolutionCache::unperturbed_values(
    const common::Range1D<Uint> elem_range) const
{
  const math::DenseConstVecView<Real> cached_values(&m_perturb_cache[elem_range.lbound()],
                                                    &m_perturb_cache[elem_range.ubound()]);
  return cached_values;
}

// ----------------------------------------------------------------------------

void SolutionCache::flush()
{
  for (typename common::DataMap<mesh::PointSetTagExt, BufferData>::iterator it = m_data_map.begin();
       it != m_data_map.end(); ++it)
  {
    (*it.data_ptr()).nb_elem_filled = 0;
  }
  m_data_index.resize(0);
}

// ----------------------------------------------------------------------------

void SolutionCache::clear()
{
  m_max_nb_blocks = 0;
  m_nb_fields     = 0;
  m_data_map.clear();
  m_perturb_cache.resize(0);
  m_local_perturb_node.resize(0);
  m_local_perturb_component.resize(0);
  m_data_index.clear();
  m_data_index.resize(0);
}

// ----------------------------------------------------------------------------

Uint SolutionCache::nb_values_in_buffer() const
{
  return m_data_index.size();
}

// ----------------------------------------------------------------------------

Uint SolutionCache::capacity() const
{
  return m_max_nb_blocks;
}

// ----------------------------------------------------------------------------

void SolutionCache::print_types() const
{
  std::cout << "=== Solution cache status ===" << std::endl;
  for (typename common::DataMap<mesh::PointSetTagExt, BufferData>::const_iterator it =
           m_data_map.cbegin();
       it != m_data_map.cend(); ++it)
  {
    const mesh::PointSetTagExt tag = it.key_value();
    std::cout << "\tType: [" << tag.std_region_tag().as_string() << "-"
              << "{key}P" << tag.key_p_order() << "-" << tag.cell_transform_id() << "-"
              << tag.local_id() << "]" << std::endl;
    std::cout << "\tNb. of elements filled: " << (*it.data_ptr()).nb_elem_filled << "/"
              << m_max_nb_blocks << std::endl;
    std::cout << "Can hold " << (*it.data_ptr()).m_values.rows() * (*it.data_ptr()).m_values.cols()
              << " values (max " << m_max_nb_blocks << " blocks, each block with " << m_nb_fields
              << " fields)" << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
