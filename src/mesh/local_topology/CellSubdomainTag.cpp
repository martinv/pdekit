#include <boost/algorithm/string.hpp>

#include "common/StringUtils.hpp"
#include "mesh/local_topology/CellSubdomainTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

CellSubdomainTag::CellSubdomainTag()
{
}

// ----------------------------------------------------------------------------

CellSubdomainTag::CellSubdomainTag(pack_type::store_type store_value) : m_tag_data(store_value)
{
}

// ----------------------------------------------------------------------------

CellSubdomainTag::CellSubdomainTag(const ElemShape elem_shape, const CellTransform adapt_op_id,
                                   const Uint pos_in_parent, const ElemShape parent_shape)
{
  m_tag_data.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(elem_shape));
  m_tag_data.set_field<1>(static_cast<std::underlying_type<CellTransform>::type>(adapt_op_id));
  m_tag_data.set_field<2>(pos_in_parent);
  m_tag_data.set_field<3>(static_cast<std::underlying_type<ElemShape>::type>(parent_shape));
}

// ----------------------------------------------------------------------------

const std::string CellSubdomainTag::as_string() const
{
  return ElemShapeInfo::name(elem_shape()) + "-" + CellTransformName::value[m_tag_data.field<1>()] +
         "-" + common::StringUtils::to_string(m_tag_data.field<2>()) + "-" +
         ElemShapeInfo::name(parent_shape());
}

// ----------------------------------------------------------------------------

const std::string CellSubdomainTag::fields_to_string(const ElemShape elem_shape,
                                                     const CellTransform adapt_op_id,
                                                     const Uint local_pos_in_parent,
                                                     const ElemShape parent_shape)
{
  return ElemShapeInfo::name(elem_shape) + "-" +
         CellTransformName::value[static_cast<std::underlying_type<CellTransform>::type>(
             adapt_op_id)] +
         "-" + common::StringUtils::to_string(local_pos_in_parent) + "-" +
         ElemShapeInfo::name(parent_shape);
}

// ----------------------------------------------------------------------------

const CellSubdomainTag CellSubdomainTag::string_to_tag(const std::string &description)
{
  std::vector<std::string> strs;
  boost::split(strs, description, boost::is_any_of("-"));

  if (strs.size() != 6)
  {
    std::cerr << "LocalIncidencePatternTag::string_to_tag: tag '" << description << "' is invalid."
              << std::endl;
  }

  CellSubdomainTag tag;

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[0])
    {
      tag.m_tag_data.set_field<0>(
          static_cast<std::underlying_type<ElemShape>::type>(ElemShapeInfo::value(i)));
      break;
    }
  }

  for (Uint i = 0; i < CellTransformName::NbInstances; ++i)
  {
    if (CellTransformName::value[i] == strs[1])
    {
      tag.m_tag_data.set_field<1>(
          static_cast<std::underlying_type<CellTransform>::type>(CellTransformValue::value[i]));
      break;
    }
  }

  tag.m_tag_data.set_field<2>(common::StringUtils::from_string<Uint>(strs[2]));

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[3])
    {
      tag.m_tag_data.set_field<3>(
          static_cast<std::underlying_type<ElemShape>::type>(ElemShapeInfo::value(i)));
      break;
    }
  }

  return tag;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const CellSubdomainTag &tag)
{
  os << tag.store_value();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
