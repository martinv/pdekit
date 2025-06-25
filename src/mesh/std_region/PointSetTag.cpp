#include <vector>

#include <boost/algorithm/string.hpp>

#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// typedef common::Tag3<ElemShape, PolyOrder, RefTopology> PointSetTag;

// ----------------------------------------------------------------------------

PointSetTag::PointSetTag() : m_tag_data(pack_type())
{
}

// ----------------------------------------------------------------------------

PointSetTag::PointSetTag(pack_type::store_type store_value) : m_tag_data(store_value)
{
}

// ----------------------------------------------------------------------------

PointSetTag::PointSetTag(const ElemShape eshape, const Uint poly_order, const PointSetID ref_topo)
{
  m_tag_data.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(eshape));
  m_tag_data.set_field<1>(poly_order);
  m_tag_data.set_field<2>(static_cast<std::underlying_type<PointSetID>::type>(ref_topo));
}

// ----------------------------------------------------------------------------

void PointSetTag::decompose_into_fields(const PointSetTag tag, ElemShape &eshape, Uint &poly_order,
                                        PointSetID &ref_topo)
{
  eshape     = tag.elem_shape();
  poly_order = tag.poly_order();
  ref_topo   = tag.ref_topology();
}

// ----------------------------------------------------------------------------

const std::string PointSetTag::as_string() const
{
  return ElemShapeInfo::name(elem_shape()) + "-" + PolyOrder::name(poly_order()) + "-" +
         PointSetInfo::name(ref_topology());
}

// ----------------------------------------------------------------------------

const std::string PointSetTag::fields_to_string(const ElemShape elem_shape, const Uint poly_order,
                                                const PointSetID ref_topology)
{
  return ElemShapeInfo::name(elem_shape) + "-" + PolyOrder::name(poly_order) + "-" +
         PointSetInfo::name(ref_topology);
}

// ----------------------------------------------------------------------------

const PointSetTag PointSetTag::string_to_tag(const std::string &description)
{
  std::vector<std::string> strs;
  boost::split(strs, description, boost::is_any_of("-"));

  PointSetTag tag;

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[0])
    {
      tag.m_tag_data.set_field<0>(
          static_cast<std::underlying_type<ElemShape>::type>(ElemShapeInfo::value(i)));
      break;
    }
  }

  for (Uint i = 0; i < PolyOrder::nb_instances(); ++i)
  {
    if (PolyOrder::name(PolyOrder::value(i)) == strs[1])
    {
      tag.m_tag_data.set_field<1>(PolyOrder::value(i));
      break;
    }
  }

  for (Uint i = 0; i < PointSetInfo::nb_instances(); ++i)
  {
    if (PointSetInfo::name(PointSetInfo::value(i)) == strs[2])
    {
      tag.m_tag_data.set_field<2>(
          static_cast<std::underlying_type<PointSetID>::type>(PointSetInfo::value(i)));
      break;
    }
  }

  return tag;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const PointSetTag &tag)
{
  os << tag.store_value();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
