#include "interpolation/PolyCache.hpp"
#include "mesh/point_set/StdPointSet.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------

PolyCache::PolyCache()
{
  m_poly_data.resize(0);
  m_cached_values.resize(0);
}

// ----------------------------------------------------------------------------

void PolyCache::allocate(const PolyCacheDescription &cache_desc)
{
  m_cached_values.reserve(cache_desc.nb_cached_entries());
  m_cached_values.resize(0);

  m_poly_data.resize(0);

  for (Uint b = 0; b < cache_desc.nb_blocks(); ++b)
  {
    create_new_poly_data(cache_desc.block_type(b));
  }
}

// ----------------------------------------------------------------------------

void PolyCache::push_back_to_buffer(const mesh::DiscreteElemKey key)
{
  for (Uint i = 0; i < m_poly_data.size(); ++i)
  {
    PolyData &pdata = m_poly_data[i];
    if (pdata.pdata_type == key)
    {
      m_cached_values.push_back(i);

      m_data_index.push_back(std::tuple<Uint, Uint>(i, pdata.nb_cached_values));
      pdata.nb_cached_values++;

      return;
    }
  }
}

// ----------------------------------------------------------------------------

void PolyCache::flush()
{
  m_cached_values.resize(0);
}

// ----------------------------------------------------------------------------

void PolyCache::clear()
{
  m_poly_data.clear();
  m_poly_data.resize(0);

  m_cached_values.resize(0);
}

// ----------------------------------------------------------------------------

Uint PolyCache::nb_values_in_buffer() const
{
  return m_cached_values.size();
}

// ----------------------------------------------------------------------------

Uint PolyCache::capacity() const
{
  return m_cached_values.capacity();
}

// ----------------------------------------------------------------------------

const interpolation::FEValues &PolyCache::fe_values(const Uint i) const
{
  return m_poly_data[m_cached_values[i]].ref_fe_values;
}

// ----------------------------------------------------------------------------

PolyCache::cellwise_view const PolyCache::cellwise_values(const Uint idx) const
{
  const common::PtrHandle<PolyData const> handle(&m_poly_data[m_cached_values[idx]]);
  return cellwise_view{handle};
}

// ----------------------------------------------------------------------------

const mesh::DiscreteElemKey PolyCache::key(const Uint idx) const
{
  const PolyData &pdata = m_poly_data[m_cached_values[idx]];
  return pdata.pdata_type;
}

// ----------------------------------------------------------------------------

Uint PolyCache::position_in_block(const Uint idx) const
{
  return std::get<POS_IN_PDATA>(m_data_index[idx]);
}

// ----------------------------------------------------------------------------

void PolyCache::print() const
{
  std::cout << "Polynomial cache:: cached polynomial types: " << std::endl;
  for (auto poly_data : m_poly_data)
  {
    const auto data_type = poly_data.pdata_type;
    std::cout << "  { [" << data_type.support() << "] ["
              << SFuncInfo::name(data_type.basis().shape_function()) << "-"
              << PolyOrder::name(data_type.basis().poly_order()) << "] [" << data_type.eval_pts()
              << "] }" << std::endl;
    std::cout << "  Number of cached values: " << poly_data.nb_cached_values << std::endl;
    // std::cout << "  Reference polynomial values: " << std::endl;
    // poly_data.ref_fe_values.print();
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

void PolyCache::create_new_poly_data(mesh::DiscreteElemKey const &block_type)
{
  PolyData poly_data;
  poly_data.pdata_type = block_type;

  const mesh::PointSetTag std_reg_tag = block_type.support().std_region_tag();

  const mesh::sf::SFTag sf_tag = block_type.basis();

  poly_data.ref_fe_values.configure(std_reg_tag, sf_tag);

  mesh::StdPointSet quad;
  quad.change_type(block_type.eval_pts().std_region_tag());
  poly_data.ref_fe_values.fill_Vandermonde(quad.get().coordinates(), quad.get().weights());

  poly_data.nb_cached_values = 0;

  m_poly_data.emplace_back(poly_data);
}

// ----------------------------------------------------------------------------

PolyCache::PolyData::PolyData()
{
  nb_cached_values = 0;
}

// ----------------------------------------------------------------------------

PolyCache::PolyData::PolyData(const PolyData &rhs)
    : pdata_type(rhs.pdata_type), nb_cached_values(rhs.nb_cached_values)
{
  interpolation::FEValues::copy(rhs.ref_fe_values, ref_fe_values);
}

// ----------------------------------------------------------------------------

PolyCache::PolyData::~PolyData()
{
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
