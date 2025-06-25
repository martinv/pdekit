#ifndef PDEKIT_Mesh_Shape_Function_Modal_Basis_Tag_hpp
#define PDEKIT_Mesh_Shape_Function_Modal_Basis_Tag_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// typedef common::Tag2<PrimeBasis, ElemShape> PrimeBasisTag;

// ----------------------------------------------------------------------------

class ModalBasisTag
{
  public:
  using pack_type = typename common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ModalBasisInfo::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value>;

  public:
  using store_type = typename pack_type::store_type;

  /// Default constructor
  ModalBasisTag();

  /// Construct from value
  ModalBasisTag(store_type store_value);

  /// Construct from prime basis and element shape
  ModalBasisTag(const ModalBasis prime_basis, const ElemShape elem_shape);

  /// Return string representation
  const std::string as_string() const;

  /// Return the raw stored value
  store_type store_value() const;

  /// Get element shape
  ModalBasis prime_basis() const;

  /// Return polynomial order
  ElemShape elem_shape() const;

  /// Comparison operator for ordering
  bool operator<(const ModalBasisTag &other) const;

  /// Generate a string description of tag from keys
  static const std::string fields_to_string(const ModalBasis prime_basis,
                                            const ElemShape elem_shape);

  /// Generate a tag based on its string description
  static const ModalBasisTag string_to_tag(const std::string &description);

  /// Print the tag value for debugging
  friend std::ostream &operator<<(std::ostream &os, const ModalBasisTag &tag);

  private:
  pack_type m_pack;
};

// ----------------------------------------------------------------------------

inline ModalBasisTag::store_type ModalBasisTag::store_value() const
{
  return m_pack.store_value();
}

// ----------------------------------------------------------------------------

inline ModalBasis ModalBasisTag::prime_basis() const
{
  return static_cast<ModalBasis>(m_pack.field<0>());
}

// ----------------------------------------------------------------------------

inline ElemShape ModalBasisTag::elem_shape() const
{
  return static_cast<ElemShape>(m_pack.field<1>());
}

// ----------------------------------------------------------------------------

inline bool ModalBasisTag::operator<(const ModalBasisTag &other) const
{
  return m_pack < other.m_pack;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
