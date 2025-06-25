#ifndef PDEKIT_Mesh_Point_Set_Tag_hpp
#define PDEKIT_Mesh_Point_Set_Tag_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"
#include <tuple>

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// typedef common::Tag3<ElemShape, PolyOrder, RefTopology> PointSetTag;

// ----------------------------------------------------------------------------

class PointSetTag
{
  private:
  using pack_type = typename common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<PolyOrder::nb_instances()>::value,
      common::MinNbBitsToStoreNumber<PointSetInfo::nb_instances()>::value + 10>;

  public:
  using store_type = typename pack_type::store_type;

  /// Default constructor
  PointSetTag();

  /// Construct from value
  PointSetTag(typename pack_type::store_type store_value);

  /// Construct from element shape, polynomial order and reference topology
  PointSetTag(const ElemShape eshape, const Uint poly_order, const PointSetID ref_topo);

  /// Decompose this tag into its fields and store them in separate integers
  static void decompose_into_fields(const PointSetTag tag, ElemShape &eshape, Uint &poly_order,
                                    PointSetID &ref_topo);

  /// Return string representation
  const std::string as_string() const;

  /// Return the raw stored value
  store_type store_value() const;

  /// Get element shape
  ElemShape elem_shape() const;

  /// Return polynomial order
  Uint poly_order() const;

  /// Return reference topology
  PointSetID ref_topology() const;

  /// Comparison operator for ordering
  bool operator<(const PointSetTag &rhs) const;

  /// Equality operator
  bool operator==(const PointSetTag &rhs) const;

  /// Negated equality operator
  bool operator!=(const PointSetTag &rhs) const;

  /// Generate a string description of tag from keys
  static const std::string fields_to_string(const ElemShape elem_shape, const Uint poly_order,
                                            const PointSetID ref_topology);

  /// Generate a tag based on its string description
  static const PointSetTag string_to_tag(const std::string &description);

  private:
  pack_type m_tag_data;
};

// ----------------------------------------------------------------------------

inline PointSetTag::store_type PointSetTag::store_value() const
{
  return m_tag_data.store_value();
}

// ----------------------------------------------------------------------------

inline ElemShape PointSetTag::elem_shape() const
{
  return static_cast<ElemShape>(m_tag_data.field<0>());
}

// ----------------------------------------------------------------------------

inline Uint PointSetTag::poly_order() const
{
  return m_tag_data.field<1>();
}

// ----------------------------------------------------------------------------

inline PointSetID PointSetTag::ref_topology() const
{
  return static_cast<PointSetID>(m_tag_data.field<2>());
}

// ----------------------------------------------------------------------------

inline bool PointSetTag::operator<(const PointSetTag &rhs) const
{
  return m_tag_data < rhs.m_tag_data;
}

// ----------------------------------------------------------------------------

inline bool PointSetTag::operator==(const PointSetTag &rhs) const
{
  return m_tag_data == rhs.m_tag_data;
}

// ----------------------------------------------------------------------------

inline bool PointSetTag::operator!=(const PointSetTag &rhs) const
{
  return m_tag_data != rhs.m_tag_data;
}

// ----------------------------------------------------------------------------

struct PointSetTagHasher
{
  std::size_t operator()(const PointSetTag &tag) const
  {
    return tag.store_value();
  }
};

// ----------------------------------------------------------------------------

struct PointSetTagEqualOp
{
  using value_type = PointSetTag;

  bool operator()(const PointSetTag &tag1, const PointSetTag &tag2) const
  {
    return (tag1.store_value() == tag2.store_value());
  }
};

// ----------------------------------------------------------------------------

/*
inline bool operator==(const PointSetTag tag_left, const PointSetTag tag_right)
{
  return ((tag_left.elem_shape() == tag_right.elem_shape()) &&
          (tag_left.poly_order() == tag_right.poly_order()) &&
          (tag_left.ref_topology() == tag_right.ref_topology()));
}
*/

// ----------------------------------------------------------------------------

class PointSetTagTupleHash
    : public std::unary_function<std::tuple<PointSetTag, PointSetTag>, std::size_t>
{
  public:
  inline std::size_t operator()(const std::tuple<PointSetTag, PointSetTag> &key) const
  {
    const std::size_t hash_val_1 = std::get<0>(key).store_value();
    const std::size_t hash_val_2 = std::get<1>(key).store_value();

    return (hash_val_1 << 10) + hash_val_2;
  }
};

// ----------------------------------------------------------------------------

// Print the tag value for debugging

std::ostream &operator<<(std::ostream &os, const PointSetTag &tag);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
