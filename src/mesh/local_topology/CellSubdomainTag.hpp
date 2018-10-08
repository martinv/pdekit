#ifndef PDEKIT_Mesh_Local_Topology_Cell_Subdomain_Tag_hpp
#define PDEKIT_Mesh_Local_Topology_Cell_Subdomain_Tag_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"
#include "mesh/CellTransform.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

// [element shape, local cell adaptation op. id, local id after adaptation,
//                                                         parent elem shape]
// - element shape is the shape of cell, facet, edge
// - cell adaptation operation determines how the 'parent' cell was split
// - local id after adaptation determines which of the subdomains obtained
//   by the splitting is considered.
// - parent elem shape is the shape of parent elements

// ----------------------------------------------------------------------------

// The element shape is the geometrical shape of the facet of the super-element
// In the picture below, the super-element would be a cube in 3D and its facet
// is a quadrilateral 6 bits for split variant on the left side of facet 5 bits
// for split variant on the right side of facet Example: there can be a split
// variant for a square facet to be divided into two smaller quads
//          If this is done for both and left side of the facet, the quads can
//          be mutually oriented as
//
//          1)
//
//          *---------------*     *---------------*
//          |               |     |               |
//          |      a1       |     |      a2       |
//          *---------------*     *---------------*
//          |               |     |               |
//          |      b1       |     |      b2       |
//          *---------------*     *---------------*
//
//          2)
//
//          *-------*-------*     *---------------*
//          |       |       |     |               |
//          |       |       |     |      a2       |
//          |   a1  |   b1  |     *---------------*
//          |       |       |     |               |
//          |       |       |     |      b2       |
//          *-------*-------*     *---------------*
//
// i.e. in the second case, the quads 'cross'
// To distinguish 1) and 2), the split variant of the left side of the facet in
// 1) and 2) has to be considered as a different case, NOT the same variant with
// rotation applied

class CellSubdomainTag
{
  private:
  typedef typename common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value, 7, 7,
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value>
      pack_type;

  public:
  typedef typename pack_type::store_type store_type;

  /// Default constructor
  CellSubdomainTag();

  /// Construct from value
  CellSubdomainTag(typename pack_type::store_type store_value);

  /// Construct from element shape, polynomial order and reference topology
  CellSubdomainTag(const ElemShape elem_shape, const CellTransform parent_adapt_op_id,
                   const Uint pos_in_parent, const ElemShape parent_shape);

  /// Return string representation
  const std::string as_string() const;

  /// Return the raw stored value
  store_type store_value() const;

  /// Get element shape
  ElemShape elem_shape() const;

  /// Get the adaptation operation id of the parent that
  /// resulted in this entity
  CellTransform parent_adapt_op_id() const;

  /// Get the index of subdomain (and hence the subdomain position)
  /// within the parent cell
  Uint local_pos_in_parent() const;

  /// Get the shape of the parent cell (whose refinement led to
  /// the creation of this cell)
  ElemShape parent_shape() const;

  /// Comparison operator for ordering
  bool operator<(const CellSubdomainTag &other) const;

  /// Generate a string description of tag from keys
  static const std::string fields_to_string(const ElemShape elem_shape,
                                            const CellTransform parent_adapt_op_id,
                                            const Uint local_pos_in_parent,
                                            const ElemShape parent_shape);

  /// Generate a tag based on its string description
  static const CellSubdomainTag string_to_tag(const std::string &description);

  private:
  pack_type m_tag_data;
};

// ----------------------------------------------------------------------------

inline CellSubdomainTag::store_type CellSubdomainTag::store_value() const
{
  return m_tag_data.store_value();
}

// ----------------------------------------------------------------------------

inline ElemShape CellSubdomainTag::elem_shape() const
{
  return static_cast<ElemShape>(m_tag_data.field<0>());
}

// ----------------------------------------------------------------------------

inline CellTransform CellSubdomainTag::parent_adapt_op_id() const
{
  return static_cast<CellTransform>(m_tag_data.field<1>());
}

// ----------------------------------------------------------------------------

inline Uint CellSubdomainTag::local_pos_in_parent() const
{
  return m_tag_data.field<2>();
}

// ----------------------------------------------------------------------------

inline ElemShape CellSubdomainTag::parent_shape() const
{
  return static_cast<ElemShape>(m_tag_data.field<3>());
}

// ----------------------------------------------------------------------------

inline bool CellSubdomainTag::operator<(const CellSubdomainTag &other) const
{
  return m_tag_data < other.m_tag_data;
}

// ----------------------------------------------------------------------------

struct FacetIncidencePatternTagHasher
{
  std::size_t operator()(const CellSubdomainTag &tag) const
  {
    return tag.store_value();
  }
};

// ----------------------------------------------------------------------------

struct FacetIncidencePatternTagEqualOp
{
  typedef CellSubdomainTag value_type;

  bool operator()(const CellSubdomainTag &tag1, const CellSubdomainTag &tag2) const
  {
    return (tag1.store_value() == tag2.store_value());
  }
};

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const CellSubdomainTag &tag);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
