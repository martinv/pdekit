#include <vector>

#include <boost/algorithm/string.hpp>

#include "mesh/shape_function//ModalBasisTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// typedef common::Tag3<ElemShape, PolyOrder, RefTopology> PrimeBasisTag;

// ----------------------------------------------------------------------------

ModalBasisTag::ModalBasisTag() : m_pack(pack_type())
{
}

// ----------------------------------------------------------------------------

ModalBasisTag::ModalBasisTag(store_type store_value) : m_pack(store_value)
{
}

// ----------------------------------------------------------------------------

ModalBasisTag::ModalBasisTag(const ModalBasis prime_basis, const ElemShape elem_shape)
{
  m_pack.set_field<0>(static_cast<std::underlying_type<ModalBasis>::type>(prime_basis));
  m_pack.set_field<1>(static_cast<std::underlying_type<ElemShape>::type>(elem_shape));
}

// ----------------------------------------------------------------------------

const std::string ModalBasisTag::as_string() const
{
  return ModalBasisInfo::name(prime_basis()) + "-" + ElemShapeInfo::name(elem_shape());
}

// ----------------------------------------------------------------------------

const std::string ModalBasisTag::fields_to_string(const ModalBasis prime_basis,
                                                  const ElemShape elem_shape)
{
  return ModalBasisInfo::name(prime_basis) + "-" + ElemShapeInfo::name(elem_shape);
}

// ----------------------------------------------------------------------------

const ModalBasisTag ModalBasisTag::string_to_tag(const std::string &description)
{
  std::vector<std::string> strs;
  boost::split(strs, description, boost::is_any_of("-"));

  ModalBasisTag tag;

  for (Uint i = 0; i < ModalBasisInfo::nb_instances(); ++i)
  {
    if (ModalBasisInfo::name(ModalBasisInfo::value(i)) == strs[0])
    {
      tag.m_pack.set_field<0>(
          static_cast<std::underlying_type<ModalBasis>::type>(ModalBasisInfo::value(i)));
      break;
    }
  }

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[1])
    {
      tag.m_pack.set_field<1>(
          static_cast<std::underlying_type<ElemShape>::type>(ElemShapeInfo::value(i)));
      break;
    }
  }

  return tag;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const ModalBasisTag &tag)
{
  os << tag.store_value();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
