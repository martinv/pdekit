#include <type_traits>

#include "common/StringUtils.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace mesh
{

// --------------------------------------------------------------------------

PointSetTagExt::PointSetTagExt()
    : m_std_reg_tag(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined))
{
  m_refinement_data.set_field<0>(P0);
  m_refinement_data.set_field<1>(static_cast<Uint>(CellTransform::NO_TRANS));
  m_refinement_data.set_field<2>(0);
}

// --------------------------------------------------------------------------

PointSetTagExt::PointSetTagExt(const PointSetTag tag, const Uint order,
                               const CellTransform cell_transform_id, const Uint loc_id)
    : m_std_reg_tag(tag)
{
  m_refinement_data.set_field<0>(order);
  m_refinement_data.set_field<1>(
      static_cast<std::underlying_type<CellTransform>::type>(cell_transform_id));
  m_refinement_data.set_field<2>(loc_id);
}

// --------------------------------------------------------------------------

PointSetTagExt::PointSetTagExt(const PointSetTag tag, const Uint order) : m_std_reg_tag(tag)
{
  m_refinement_data.set_field<0>(order);
  m_refinement_data.set_field<1>(
      static_cast<std::underlying_type<CellTransform>::type>(CellTransform::NO_TRANS));
  m_refinement_data.set_field<2>(0u);
}

// --------------------------------------------------------------------------

PointSetTagExt::~PointSetTagExt()
{
}

// --------------------------------------------------------------------------

const std::string PointSetTagExt::as_string() const
{
  return m_std_reg_tag.as_string() + "-{key}P" + common::StringUtils::to_string(key_p_order()) +
         "-" + common::StringUtils::to_string(cell_transform_id()) + "-" +
         common::StringUtils::to_string(local_id());
}

// --------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, PointSetTagExt const &key)
{
  os << "[" << key.std_region_tag().as_string() << "-{key}P" << key.key_p_order() << "-"
     << key.cell_transform_id() << "-" << key.local_id() << "]";
  return os;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const PointSetTagExtPair &tag_pair)
{
  os << "[" << tag_pair.m_t1 << "," << tag_pair.m_t2 << "]";
  return os;
}

// --------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
