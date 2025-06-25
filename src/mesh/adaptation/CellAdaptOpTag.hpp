#ifndef PDEKIT_Mesh_Adaptation_Cell_Adapt_Op_Tag_hpp
#define PDEKIT_Mesh_Adaptation_Cell_Adapt_Op_Tag_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"
#include "mesh/CellTransform.hpp"

namespace pdekit
{

namespace mesh
{

namespace adapt
{

// ----------------------------------------------------------------------------

// [ElemShape eshape, Uint adapt_op_id]

// ----------------------------------------------------------------------------

class CellAdaptOpTag
{
  private:
  typedef typename common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value, 11>
      pack_type;

  public:
  typedef typename pack_type::store_type store_type;

  /// Default constructor
  CellAdaptOpTag();

  /// Construct from value
  CellAdaptOpTag(store_type store_value);

  /// Construct from element shape, polynomial order and reference topology
  CellAdaptOpTag(const ElemShape eshape, const CellTransform adapt_op_id);

  /// Decompose this tag into its fields and store them in separate integers
  static void decompose_into_fields(const CellAdaptOpTag tag, ElemShape &eshape,
                                    CellTransform &adapt_op_id);

  /// Return string representation
  const std::string as_string() const;

  /// Set element shape
  // void set_elem_shape(const ElemShape eshape);

  /// Set type of adaptation operation
  // void set_adapt_op_id(const CellTransform adapt_op_id);

  ElemShape elem_shape() const;

  /// Return split variant on the left side of the facet
  CellTransform adapt_op_id() const;

  /// Return the raw stored value
  store_type store_value() const;

  /// Equality operator
  bool operator==(const CellAdaptOpTag &rhs) const;

  /// Comparison operator for ordering
  bool operator<(const CellAdaptOpTag &other) const;

  /// Generate a string description of tag from keys
  static const std::string fields_to_string(const ElemShape elem_shape,
                                            const CellTransform adapt_op_id);

  /// Generate a tag based on its string description
  static const CellAdaptOpTag string_to_tag(const std::string &description);

  /// Print the tag value for debugging
  friend std::ostream &operator<<(std::ostream &os, const CellAdaptOpTag &tag);

  private:
  pack_type m_pack;
};

// ----------------------------------------------------------------------------

inline ElemShape CellAdaptOpTag::elem_shape() const
{
  return static_cast<ElemShape>(m_pack.field<0>());
}

// ----------------------------------------------------------------------------

inline CellTransform CellAdaptOpTag::adapt_op_id() const
{
  return static_cast<CellTransform>(m_pack.field<1>());
}

// ----------------------------------------------------------------------------

inline CellAdaptOpTag::store_type CellAdaptOpTag::store_value() const
{
  return m_pack.store_value();
}

// ----------------------------------------------------------------------------

inline bool CellAdaptOpTag::operator==(const CellAdaptOpTag &other) const
{
  return m_pack == other.m_pack;
}

// ----------------------------------------------------------------------------

inline bool CellAdaptOpTag::operator<(const CellAdaptOpTag &other) const
{
  return m_pack < other.m_pack;
}

// ----------------------------------------------------------------------------

struct CellSplitStrategyTagHasher
{
  std::size_t operator()(const CellAdaptOpTag &tag) const
  {
    return tag.store_value();
  }
};

// ----------------------------------------------------------------------------

struct CellSplitStrategyTagEqualOp
{
  typedef CellAdaptOpTag value_type;

  bool operator()(const CellAdaptOpTag &tag1, const CellAdaptOpTag &tag2) const
  {
    return (tag1.store_value() == tag2.store_value());
  }
};

// ----------------------------------------------------------------------------

} // namespace adapt

} // namespace mesh

} // namespace pdekit

#endif
