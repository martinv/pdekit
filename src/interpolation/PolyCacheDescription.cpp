#include "interpolation/PolyCacheDescription.hpp"

namespace pdekit
{

namespace interpolation
{

// ----------------------------------------------------------------------------
// Methods of PolyCacheDescription
// ----------------------------------------------------------------------------

PolyCacheDescription::PolyCacheDescription()
{
  m_blocks.resize(0);
}

// ----------------------------------------------------------------------------

PolyCacheDescription::PolyCacheDescription(const PolyCacheDescription &rhs)
{
  m_blocks.resize(rhs.m_blocks.size());
  std::copy(rhs.m_blocks.begin(), rhs.m_blocks.end(), m_blocks.begin());
}

// ----------------------------------------------------------------------------

PolyCacheDescription::~PolyCacheDescription()
{
}

// ----------------------------------------------------------------------------

const PolyCacheDescription &PolyCacheDescription::operator=(const PolyCacheDescription &rhs)
{
  m_blocks.resize(rhs.m_blocks.size());
  std::copy(rhs.m_blocks.begin(), rhs.m_blocks.end(), m_blocks.begin());
  return *this;
}

// ----------------------------------------------------------------------------

void PolyCacheDescription::add_to_block_count(const mesh::DiscreteElemKey &type,
                                              const Uint nb_values)
{
  for (auto bdata : m_blocks)
  {
    if (std::get<0>(bdata) == type)
    {
      std::get<1>(bdata) += nb_values;
      return;
    }
  }

  std::tuple<mesh::DiscreteElemKey, Uint> new_bdata(type, nb_values);
  // std::get<0>(new_bdata) = type;
  // std::get<1>(new_bdata) = 1;
  m_blocks.push_back(new_bdata);
}

// ----------------------------------------------------------------------------

void PolyCacheDescription::add_to_block_count(const mesh::PointSetTagExt &std_reg_key,
                                              const mesh::sf::SFTag sf_type,
                                              const mesh::PointSetTagExt quad_tag,
                                              const Uint nb_values)
{
  const mesh::DiscreteElemKey type(std_reg_key, sf_type, quad_tag);
  for (auto bdata : m_blocks)
  {
    if (std::get<0>(bdata) == type)
    {
      std::get<1>(bdata) += nb_values;
      return;
    }
  }

  std::tuple<mesh::DiscreteElemKey, Uint> new_bdata(type, nb_values);
  // std::get<0>(new_bdata) = type;
  // std::get<1>(new_bdata) = 1;
  m_blocks.push_back(new_bdata);
}

// ----------------------------------------------------------------------------

Uint PolyCacheDescription::block_size(const mesh::DiscreteElemKey &type) const
{
  for (auto bdata : m_blocks)
  {
    if (std::get<0>(bdata) == type)
    {
      return std::get<1>(bdata);
    }
  }
  return 0;
}

// ----------------------------------------------------------------------------

std::vector<std::tuple<mesh::DiscreteElemKey, Uint>> PolyCacheDescription::all_block_sizes() const
{
  std::vector<std::tuple<mesh::DiscreteElemKey, Uint>> blocks(m_blocks);
  return blocks;
}

// ----------------------------------------------------------------------------

const mesh::DiscreteElemKey &PolyCacheDescription::block_type(const Uint block_idx) const
{
  return std::get<0>(m_blocks[block_idx]);
}

// ----------------------------------------------------------------------------

Uint PolyCacheDescription::nb_blocks() const
{
  return m_blocks.size();
}

// ----------------------------------------------------------------------------

Uint PolyCacheDescription::nb_cached_entries() const
{
  Uint result = 0;
  for (auto bdata : m_blocks)
  {
    result += std::get<1>(bdata);
  }
  return result;
}

// ----------------------------------------------------------------------------

void PolyCacheDescription::clear()
{
  m_blocks.resize(0);
}

// ----------------------------------------------------------------------------

void PolyCacheDescription::print() const
{
  std::cout << "Polynomial cache description contents:" << std::endl;
  for (auto poly_data : m_blocks)
  {
    const mesh::DiscreteElemKey &block_type = std::get<0>(poly_data);

    std::cout << "  { [" << block_type.support() << "] ["
              << SFuncInfo::name(block_type.basis().shape_function()) << "-"
              << PolyOrder::name(block_type.basis().poly_order()) << "] [" << block_type.eval_pts()
              << "] }" << std::endl;
    std::cout << "  Number of cached values: " << std::get<1>(poly_data) << std::endl;
    // std::cout << "  Reference polynomial values: " << std::endl;
    // poly_data.ref_fe_values.print();
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

} // namespace interpolation

} // namespace pdekit
