#include <vector>

#include <boost/algorithm/string.hpp>

#include "mesh/shape_function/SFTag.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

SFTag::SFTag() : m_pack()
{
}

// ----------------------------------------------------------------------------

SFTag::SFTag(pack_type::store_type store_value) : m_pack(store_value)
{
}

// ----------------------------------------------------------------------------

SFTag::SFTag(const ElemShape eshape, const SFunc sfunc, const Uint poly_order,
             const ModalBasis prime_basis)
{
  m_pack.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(eshape));
  m_pack.set_field<1>(static_cast<std::underlying_type<SFunc>::type>(sfunc));
  m_pack.set_field<2>(poly_order);
  m_pack.set_field<3>(static_cast<std::underlying_type<ModalBasis>::type>(prime_basis));
}

// ----------------------------------------------------------------------------

void SFTag::decompose_into_fields(const SFTag tag, ElemShape &eshape, SFunc &sfunc,
                                  Uint &poly_order, ModalBasis &prime_basis)
{
  /*
  eshape = static_cast<ElemShape>(tag.field<0>());
  sfunc = static_cast<SFunc>(tag.field<1>());
  ref_topo = static_cast<PointSetID>(tag.field<2>());
  poly_order = tag.field<3>();
  prime_basis = static_cast<PrimeBasis>(tag.field<4>());
  */

  eshape      = tag.elem_shape();
  sfunc       = tag.shape_function();
  poly_order  = tag.poly_order();
  prime_basis = tag.prime_basis();
}

// ----------------------------------------------------------------------------

const std::string SFTag::as_string() const
{
  return ElemShapeInfo::name(elem_shape()) + "-" + SFuncInfo::name(shape_function()) + "-" +
         PolyOrder::name(poly_order()) + "-" + ModalBasisInfo::name(prime_basis());
}

// ----------------------------------------------------------------------------

void SFTag::set_elem_shape(const ElemShape eshape)
{
  m_pack.set_field<0>(static_cast<std::underlying_type<ElemShape>::type>(eshape));
}

// ----------------------------------------------------------------------------

void SFTag::set_shape_function(const SFunc sfunc)
{
  m_pack.set_field<1>(static_cast<std::underlying_type<SFunc>::type>(sfunc));
}

// ----------------------------------------------------------------------------

void SFTag::set_poly_order(const Uint poly_order)
{
  m_pack.set_field<2>(poly_order);
}

// ----------------------------------------------------------------------------

void SFTag::set_prime_basis(const ModalBasis prime_basis)
{
  m_pack.set_field<3>(static_cast<std::underlying_type<ModalBasis>::type>(prime_basis));
}

// ----------------------------------------------------------------------------

const std::string SFTag::fields_to_string(const ElemShape elem_shape, const SFunc shape_function,
                                          const Uint poly_order, const ModalBasis prime_basis)
{
  return ElemShapeInfo::name(elem_shape) + "-" + SFuncInfo::name(shape_function) + "-" +
         PolyOrder::name(poly_order) + "-" + ModalBasisInfo::name(prime_basis);
}

// ----------------------------------------------------------------------------

const SFTag SFTag::string_to_tag(const std::string &description)
{
  std::vector<std::string> strs;
  boost::split(strs, description, boost::is_any_of("-"));

  SFTag tag;

  for (Uint i = 0; i < ElemShapeInfo::nb_instances(); ++i)
  {
    if (ElemShapeInfo::name(ElemShapeInfo::value(i)) == strs[0])
    {
      tag.set_elem_shape(ElemShapeInfo::value(i));
      break;
    }
  }

  for (Uint i = 0; i < SFuncInfo::nb_instances(); ++i)
  {
    if (SFuncInfo::name(SFuncInfo::value(i)) == strs[1])
    {
      tag.set_shape_function(SFuncInfo::value(i));
      break;
    }
  }

  for (Uint i = 0; i < PolyOrder::nb_instances(); ++i)
  {
    if (PolyOrder::name(PolyOrder::value(i)) == strs[2])
    {
      tag.set_poly_order(PolyOrder::value(i));
      break;
    }
  }

  for (Uint i = 0; i < ModalBasisInfo::nb_instances(); ++i)
  {
    if (ModalBasisInfo::name(ModalBasisInfo::value(i)) == strs[3])
    {
      tag.set_prime_basis(ModalBasisInfo::value(i));
      break;
    }
  }

  return tag;
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const SFTag &tag)
{
  os << tag.m_pack.store_value();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit
