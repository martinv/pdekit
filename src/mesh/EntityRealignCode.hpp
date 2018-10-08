#ifndef PDEKIT_Mesh_Entity_Realign_Code_hpp
#define PDEKIT_Mesh_Entity_Realign_Code_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"
#include "mesh/CellTransform.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// This class holds information about the permutation (vertex
// reordering) of entity: how many rotations, entity flips
// and potentially some other operations are needed to put
// the entity into the new ('permuted') position
// ----------------------------------------------------------------------------

class EntityRealignCode
{
  public:
  /// TYPEDEFS
  using pack_type = common::SmallIntegerPack<
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value, 7, 7,
      common::MinNbBitsToStoreNumber<ElemShapeInfo::nb_instances()>::value, 1, 5>;

  using store_type = typename pack_type::store_type;

  /// METHODS

  /// Return identity permutation code
  static const EntityRealignCode identity(const ElemShape eshape);

  /// Return a permutation code representing a single flip
  /// We cannot return a const object, because it will be first constructed
  /// and then modified in the static method
  static const EntityRealignCode single_flip(const ElemShape eshape);

  /// Return a permutation code representing a single rotation
  /// We cannot return a const object, because it will be first constructed
  /// and then modified in the static method
  static const EntityRealignCode single_rotation(const ElemShape eshape);

  /// Default constructor
  EntityRealignCode();

  /// Construct from given permutation value
  EntityRealignCode(pack_type::store_type store_value);

  /// Construct from all fields needed to define the object
  EntityRealignCode(const ElemShape elem_shape, const CellTransform adapt_op_id,
                    const Uint pos_in_parent, const ElemShape parent_shape, const Uint nb_flips,
                    const Uint nb_rotations);

  /// Set everything to zero
  void reset();

  /// Return the raw stored value
  store_type store_value() const;

  /// Return true if the permutation is identity
  bool is_identity(const ElemShape eshape) const;

  /// Get element shape
  ElemShape elem_shape() const;

  /// Get the adaptation operation id
  CellTransform adapt_op_id() const;

  /// Get the index of subdomain (and hence the subdomain position)
  /// within the parent cell
  Uint local_pos_in_parent() const;

  /// Get the shape of the parent cell (whose refinement led to
  /// the creation of this cell)
  ElemShape parent_shape() const;

  /// Set the number of flips in this permutation
  void set_nb_flips(const Uint nb_flips);

  /// Set the number of roatations in this permutation
  void set_nb_rotations(const Uint nb_rotations);

  /// Increment the number of flips by one
  void add_flip();

  /// Increment the number of rotations by one
  void add_rotation();

  /// Get the number of flips
  Uint nb_flips() const;

  /// Get the number of rotations
  Uint nb_rotations() const;

  /// Comparison operator for ordering
  bool operator<(const EntityRealignCode &other) const;

  // Print full human-readable description of the code
  const std::string as_string() const;

  /// Print to output stream
  friend std::ostream &operator<<(std::ostream &os, EntityRealignCode const &code);

  private:
  pack_type m_tag_data;
};

// ----------------------------------------------------------------------------

inline EntityRealignCode::store_type EntityRealignCode::store_value() const
{
  return m_tag_data.store_value();
}

// ----------------------------------------------------------------------------

inline ElemShape EntityRealignCode::elem_shape() const
{
  return static_cast<ElemShape>(m_tag_data.field<0>());
}

// ----------------------------------------------------------------------------

inline CellTransform EntityRealignCode::adapt_op_id() const
{
  return static_cast<CellTransform>(m_tag_data.field<1>());
}

// ----------------------------------------------------------------------------

inline Uint EntityRealignCode::local_pos_in_parent() const
{
  return m_tag_data.field<2>();
}

// ----------------------------------------------------------------------------

inline ElemShape EntityRealignCode::parent_shape() const
{
  return static_cast<ElemShape>(m_tag_data.field<3>());
}

// ----------------------------------------------------------------------------

inline Uint EntityRealignCode::nb_flips() const
{
  return m_tag_data.field<4>();
}

// ----------------------------------------------------------------------------

inline Uint EntityRealignCode::nb_rotations() const
{
  return m_tag_data.field<5>();
}

// ----------------------------------------------------------------------------

inline bool EntityRealignCode::operator<(const EntityRealignCode &other) const
{
  return m_tag_data.store_value() < other.m_tag_data.store_value();
}

// ----------------------------------------------------------------------------

inline bool operator==(const EntityRealignCode tag_left, const EntityRealignCode tag_right)
{
  return (tag_left.store_value() == tag_right.store_value());
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
