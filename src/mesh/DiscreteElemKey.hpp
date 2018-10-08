#ifndef PDEKIT_Mesh_Discrete_Elem_Key_hpp
#define PDEKIT_Mesh_Discrete_Elem_Key_hpp

#include <vector>

#include "mesh/shape_function/SFTag.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class DiscreteElemKey
{
  public:
  DiscreteElemKey() = default;

  DiscreteElemKey(const PointSetTagExt &support, const sf::SFTag &basis,
                  const PointSetTagExt &eval_pts);

  DiscreteElemKey(const DiscreteElemKey &other) = default;

  ~DiscreteElemKey() = default;

  DiscreteElemKey &operator=(const DiscreteElemKey &rhs) = default;

  PointSetTagExt support() const;

  sf::SFTag basis() const;

  PointSetTagExt eval_pts() const;

  private:
  PointSetTagExt m_support;

  sf::SFTag m_basis;

  PointSetTagExt m_eval_pts;
};

// ----------------------------------------------------------------------------

inline PointSetTagExt DiscreteElemKey::support() const
{
  return m_support;
}

// ----------------------------------------------------------------------------

inline sf::SFTag DiscreteElemKey::basis() const
{
  return m_basis;
}

// ----------------------------------------------------------------------------

inline PointSetTagExt DiscreteElemKey::eval_pts() const
{
  return m_eval_pts;
}

// ----------------------------------------------------------------------------

class DiscreteElemKeyHash : public std::unary_function<DiscreteElemKey, std::size_t>
{
  public:
  inline std::size_t operator()(const DiscreteElemKey &key) const
  {
    const mesh::PointSetTagExt key_part0 = key.support();
    const mesh::sf::SFTag key_part1      = key.basis();
    const mesh::PointSetTagExt key_part2 = key.eval_pts();

    const std::size_t hash_val_0 =
        key_part0.std_region_tag().store_value() ^ key_part0.local_id() ^ key_part0.key_p_order();
    const std::size_t hash_val_1 = key_part1.store_value();
    const std::size_t hash_val_2 =
        key_part2.std_region_tag().store_value() ^ key_part2.local_id() ^ key_part2.key_p_order();

    return (hash_val_0 << 10) + (hash_val_1 << 5) + hash_val_2;
  }
};

// ----------------------------------------------------------------------------

inline bool operator==(const DiscreteElemKey &key_L, const DiscreteElemKey &key_R)
{
  return ((key_L.support() == key_R.support()) && (key_L.basis() == key_R.basis()) &&
          (key_L.eval_pts() == key_R.eval_pts()));
}

// ----------------------------------------------------------------------------

inline bool operator!=(const DiscreteElemKey &key_L, const DiscreteElemKey &key_R)
{
  return !operator==(key_L, key_R);
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const DiscreteElemKey &key);

// ----------------------------------------------------------------------------

void add_unique_discr_elem_key(std::vector<DiscreteElemKey> &keys,
                               const DiscreteElemKey &candidate_key);

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
