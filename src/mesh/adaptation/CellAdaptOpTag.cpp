#include <boost/algorithm/string.hpp>

#include "common/StringUtils.hpp"
#include "mesh/adaptation/CellAdaptOpTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

// ----------------------------------------------------------------------------

// typedef common::Tag3<ElemShape, PolyOrder, RefTopology> CellSplitStrategyTag;

// ----------------------------------------------------------------------------

CellAdaptOpTag::CellAdaptOpTag() : m_pack(pack_type())
{
}

// ----------------------------------------------------------------------------

CellAdaptOpTag::CellAdaptOpTag(store_type store_value) : m_pack(store_value)
{
}

// ----------------------------------------------------------------------------

CellAdaptOpTag::CellAdaptOpTag(const ElemShape eshape, const CellTransform adapt_op_id)
{
  m_pack.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(eshape));
  m_pack.set_field<1>(static_cast<std::underlying_type<CellTransform>::type>(adapt_op_id));
}

// ----------------------------------------------------------------------------

void CellAdaptOpTag::decompose_into_fields(const CellAdaptOpTag tag, ElemShape &eshape,
                                           CellTransform &adapt_op_id)
{
  eshape      = tag.elem_shape();
  adapt_op_id = tag.adapt_op_id();
}

// ----------------------------------------------------------------------------

const std::string CellAdaptOpTag::as_string() const
{
  return ElemShapeInfo::name(elem_shape()) + "-" + CellTransformName::value[m_pack.field<1>()];
}

// ----------------------------------------------------------------------------

/*
void CellAdaptOpTag::set_elem_shape(const ElemShape eshape)
{
  m_pack.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(eshape));
}

// ----------------------------------------------------------------------------

void CellAdaptOpTag::set_adapt_op_id(const CellTransform adapt_op_id)
{
  m_pack.set_field<1>(static_cast<std::underlying_type<CellTransform>::type>(adapt_op_id));
}
*/

// ----------------------------------------------------------------------------

const std::string CellAdaptOpTag::fields_to_string(const ElemShape elem_shape,
                                                   const CellTransform adapt_op_id)
{
  return ElemShapeInfo::name(elem_shape) + "-" +
         CellTransformName::value[static_cast<std::underlying_type<CellTransform>::type>(
             adapt_op_id)];
}

// ----------------------------------------------------------------------------

const CellAdaptOpTag CellAdaptOpTag::string_to_tag(const std::string &description)
{
  std::vector<std::string> strs;
  boost::split(strs, description, boost::is_any_of("-"));

  CellAdaptOpTag tag;

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[0])
    {
      tag.m_pack.set_field<0>(
          static_cast<std::underlying_type<ElemShape>::type>(ElemShapeInfo::value(i)));
      break;
    }
  }

  for (Uint i = 0; i < CellTransformName::NbInstances; ++i)
  {
    if (CellTransformName::value[i] == strs[1])
    {
      tag.m_pack.set_field<1>(
          static_cast<std::underlying_type<CellTransform>::type>(CellTransformValue::value[i]));
      break;
    }
  }

  return tag;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const CellAdaptOpTag &tag)
{
  os << tag.store_value();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit
