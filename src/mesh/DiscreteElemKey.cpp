#include "mesh/DiscreteElemKey.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

DiscreteElemKey::DiscreteElemKey(const PointSetTagExt &support, const sf::SFTag &basis,
                                 const PointSetTagExt &eval_pts)
    : m_support(support), m_basis(basis), m_eval_pts(eval_pts)
{
}

// ----------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &os, const DiscreteElemKey &key)
{
  os << "{" << key.support() << "-" << key.basis() << "-" << key.eval_pts() << "}";
  return os;
}

// ----------------------------------------------------------------------------

void add_unique_discr_elem_key(std::vector<DiscreteElemKey> &keys,
                               const DiscreteElemKey &candidate_key)
{
  bool candidate_found = false;
  for (const auto &key : keys)
  {
    if (candidate_key == key)
    {
      candidate_found = true;
      break;
    }
  }

  if (!candidate_found)
  {
    keys.push_back(candidate_key);
  }
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
