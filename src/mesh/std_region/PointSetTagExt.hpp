#ifndef PDEKIT_Mesh_Point_Set_Tag_Ext_hpp
#define PDEKIT_Mesh_Point_Set_Tag_Ext_hpp

#include "common/Constants.hpp"
#include "common/SmallIntegerPack.hpp"
#include "mesh/CellTransform.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

namespace detail
{

class PointSetTagExtHasher;
class PointSetTagExtPairHasher;

} // namespace detail

// ----------------------------------------------------------------------------

class PointSetTagExt
{
  public:
  using hash_type = detail::PointSetTagExtHasher;

  /// Default constructor
  PointSetTagExt();

  /// Constructor with value initialization
  /// @param tag               ... standard region tag - describes a reference
  /// region
  /// @param order             ... polynomial order - enables to distinguish
  /// between
  ///                              identical standard region tags which are
  ///                              associated with different polynomial orders
  /// @param cell_transform_id ... type of transformation (refinement) of
  /// parent standard region
  ///                              that generated the current standard region
  /// @param loc_id            ... local index within parent region - defines
  /// which child
  ///                              within parent standard region is this
  PointSetTagExt(const PointSetTag tag, const Uint order, const CellTransform cell_transform_id,
                 const Uint loc_id);

  /// Constructor with value initialization, suppose that
  /// only order is given and cell_transform_id = CellTransform::DO_NOTHING
  /// and loc_id = 0
  /// @param tag    ... standard region tag - describes a reference region
  /// @param order  ... polynomial order - enables to distinguish between
  ///                   identical standard region tags which are associated
  ///                   with different polynomial orders
  PointSetTagExt(const PointSetTag tag, const Uint order);

  /// Destructor
  ~PointSetTagExt();

  /// Return tag of standard region
  const PointSetTag std_region_tag() const;

  /// Return polynomial order. This polynomial order
  /// is not associated with the standard region tag in this key, but is an
  /// extra variable so that multiple DefaultStdRegMapKeys differing only by
  /// polynomial order can be used
  Uint key_p_order() const;
  /// Return the adaptation operation id

  CellTransform cell_transform_id() const;

  /// Return local id
  Uint local_id() const;

  /// Return string representation
  const std::string as_string() const;

  private:
  /// Tag for element type
  PointSetTag m_std_reg_tag;

  /// Additional data needed when dealing with refined mesh:
  /// <type of adaptation operation, polynomial order of the key, local id in
  /// parent face>
  ///
  /// The polynomial order here is not necessarily the same
  /// as the order of the std. region stored by the tag
  /// It can be quadrature order for example
  common::SmallIntegerPack<common::MinNbBitsToStoreNumber<PolyOrder::nb_instances()>::value,
                           common::MinNbBitsToStoreNumber<CellTransformValue::NbInstances>::value,
                           7>
      m_refinement_data;
};

// --------------------------------------------------------------------------

inline const PointSetTag PointSetTagExt::std_region_tag() const
{
  return m_std_reg_tag;
}

// --------------------------------------------------------------------------

inline Uint PointSetTagExt::key_p_order() const
{
  return m_refinement_data.field<0>();
}

// ----------------------------------------------------------------------------

inline CellTransform PointSetTagExt::cell_transform_id() const
{
  return static_cast<CellTransform>(m_refinement_data.field<1>());
}

// ----------------------------------------------------------------------------

inline Uint PointSetTagExt::local_id() const
{
  return m_refinement_data.field<2>();
}

// ----------------------------------------------------------------------------

inline bool operator==(const PointSetTagExt &lhs, const PointSetTagExt &rhs)
{
  return (lhs.std_region_tag() == rhs.std_region_tag()) && (lhs.local_id() == rhs.local_id()) &&
         (lhs.key_p_order() == rhs.key_p_order()) &&
         (lhs.cell_transform_id() == rhs.cell_transform_id());
}

// ----------------------------------------------------------------------------

inline bool operator<(const PointSetTagExt &lhs, const PointSetTagExt &rhs)
{
  if (lhs.std_region_tag() != rhs.std_region_tag())
  {
    return lhs.std_region_tag() < rhs.std_region_tag();
  }
  if (lhs.key_p_order() != rhs.key_p_order())
  {
    return lhs.key_p_order() < rhs.key_p_order();
  }
  if (lhs.cell_transform_id() != rhs.cell_transform_id())
  {
    return lhs.cell_transform_id() < rhs.cell_transform_id();
  }

  return lhs.local_id() < rhs.local_id();
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, PointSetTagExt const &key);

// ----------------------------------------------------------------------------

struct PointSetTagExtPair
{
  PointSetTagExtPair() : m_t1(), m_t2()
  {
  }

  PointSetTagExtPair(const mesh::PointSetTagExt t1, const mesh::PointSetTagExt t2)
      : m_t1(t1), m_t2(t2)
  {
  }

  using hash_type = detail::PointSetTagExtPairHasher;

  mesh::PointSetTagExt m_t1;
  mesh::PointSetTagExt m_t2;
};

// ----------------------------------------------------------------------------

inline bool operator==(const PointSetTagExtPair &lhs, const PointSetTagExtPair &rhs)
{
  return (lhs.m_t1 == rhs.m_t1) && (lhs.m_t2 == rhs.m_t2);
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const PointSetTagExtPair &tag_pair);

// ----------------------------------------------------------------------------

namespace detail
{

// Hash function for use with std::unordered_map
class PointSetTagExtHasher
{
  public:
  inline std::size_t operator()(const PointSetTagExt &key) const
  {
    return key.std_region_tag().store_value() ^ key.local_id() ^ key.key_p_order() ^
           static_cast<Uint>(key.cell_transform_id());
  }
};

// ----------------------------------------------------------------------------

class PointSetTagExtPairHasher
{
  public:
  inline std::size_t operator()(const PointSetTagExtPair &key) const
  {
    const std::size_t hash_val_1 = key.m_t1.std_region_tag().store_value() ^ key.m_t1.local_id() ^
                                   key.m_t1.key_p_order() ^
                                   static_cast<Uint>(key.m_t1.cell_transform_id());

    const std::size_t hash_val_2 = key.m_t2.std_region_tag().store_value() ^ key.m_t2.local_id() ^
                                   key.m_t2.key_p_order() ^
                                   static_cast<Uint>(key.m_t2.cell_transform_id());

    return (hash_val_1 << 10) + hash_val_2;
  }
};

// ----------------------------------------------------------------------------

class PointSetTagExtTupleHash
    : public std::unary_function<std::tuple<PointSetTagExt, PointSetTagExt>, std::size_t>
{
  public:
  inline std::size_t operator()(const std::tuple<PointSetTagExt, PointSetTagExt> &key) const
  {
    const PointSetTagExt pste1 = std::get<0>(key);
    const PointSetTagExt pste2 = std::get<1>(key);

    const std::size_t hash_val_1 =
        pste1.std_region_tag().store_value() ^ pste1.local_id() ^ pste1.key_p_order();
    const std::size_t hash_val_2 =
        pste2.std_region_tag().store_value() ^ pste2.local_id() ^ pste2.key_p_order();

    return (hash_val_1 << 10) + hash_val_2;
  }
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
