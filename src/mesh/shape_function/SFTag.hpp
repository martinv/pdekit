#ifndef PDEKIT_Interpolation_Shape_Function_Tag_hpp
#define PDEKIT_Interpolation_Shape_Function_Tag_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"

namespace pdekit
{

namespace mesh
{

namespace sf
{

// ----------------------------------------------------------------------------

// typedef common::Tag5<ElemShape, SFunc, PolyOrder, ModalBasis> SFTag;

// ----------------------------------------------------------------------------

class SFTag
{
  private:
  using pack_type = typename common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<SFuncInfo::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<PolyOrder::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<ModalBasisInfo::nb_instances()>::value + 7u>;

  public:
  using store_type = typename pack_type::store_type;

  /// Default constructor
  SFTag();

  /// Construct from value
  SFTag(typename pack_type::store_type store_value);

  /// Construct from element shape, polynomial order and reference topology
  SFTag(const ElemShape eshape, const SFunc sfunc, const Uint poly_order,
        const ModalBasis prime_basis);

  /// Decompose this tag into its fields and store them in separate integers
  static void decompose_into_fields(const SFTag tag, ElemShape &eshape, SFunc &sfunc,
                                    Uint &poly_order, ModalBasis &prime_basis);

  /// Return string representation
  const std::string as_string() const;

  /// Set element shape
  void set_elem_shape(const ElemShape eshape);

  /// Set shape function
  void set_shape_function(const SFunc sfunc);

  /// Set polynomial order
  void set_poly_order(const Uint poly_order);

  /// Set prime basis
  void set_prime_basis(const ModalBasis prime_basis);

  /// Return the raw stored value
  store_type store_value() const;

  /// Get element shape
  ElemShape elem_shape() const;

  /// Return shape function type
  SFunc shape_function() const;

  /// Return polynomial order
  Uint poly_order() const;

  /// Return prime basis type
  ModalBasis prime_basis() const;

  /// Equality operator
  bool operator==(const SFTag &rhs) const;

  /// Comparison operator for ordering
  bool operator<(const SFTag &other) const;

  /// Generate a string description of tag from keys
  static const std::string fields_to_string(const ElemShape elem_shape, const SFunc shape_function,
                                            const Uint poly_order, const ModalBasis prime_basis);

  /// Generate a tag based on its string description
  static const SFTag string_to_tag(const std::string &description);

  /// Print the tag value for debugging
  friend std::ostream &operator<<(std::ostream &os, const SFTag &tag);

  public:
  pack_type m_pack;
};

// ----------------------------------------------------------------------------

inline SFTag::store_type SFTag::store_value() const
{
  return m_pack.store_value();
}

// ----------------------------------------------------------------------------

inline ElemShape SFTag::elem_shape() const
{
  return static_cast<ElemShape>(m_pack.field<0>());
}

// ----------------------------------------------------------------------------

inline SFunc SFTag::shape_function() const
{
  return static_cast<SFunc>(m_pack.field<1>());
}

// ----------------------------------------------------------------------------

inline Uint SFTag::poly_order() const
{
  return m_pack.field<2>();
}

// ----------------------------------------------------------------------------

inline ModalBasis SFTag::prime_basis() const
{
  return static_cast<ModalBasis>(m_pack.field<3>());
}

// ----------------------------------------------------------------------------

inline bool SFTag::operator==(const SFTag &other) const
{
  return m_pack == other.m_pack;
}

// ----------------------------------------------------------------------------

inline bool SFTag::operator<(const SFTag &other) const
{
  return m_pack < other.m_pack;
}

// ----------------------------------------------------------------------------

} // namespace sf

} // namespace mesh

} // namespace pdekit

#endif
