#ifndef PDEKIT_Mesh_Entity_Dof_Realign_hpp
#define PDEKIT_Mesh_Entity_Dof_Realign_hpp

#include <iostream>
#include <vector>

#include "common/Flyweight.hpp"
#include "mesh/EntityRealignCode.hpp"
#include "mesh/std_region/PointSetTag.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

class EntityDofRealignInstance
{
  public:
  /// STATIC VARIABLES

  static const std::pair<PointSetTag, EntityRealignCode> undefined;

  /// Static method to fill a permutation
  static void construct(const std::pair<PointSetTag, EntityRealignCode> &key,
                        EntityDofRealignInstance &entity_permutation);

  /// Default constructor
  EntityDofRealignInstance();

  /// Constructor permutation given interpolation point set id and permutation
  /// code
  EntityDofRealignInstance(const std::pair<PointSetTag, EntityRealignCode> &key);

  /// Copy constructor
  EntityDofRealignInstance(const EntityDofRealignInstance &other_permutation);

  /// Assignment operator
  EntityDofRealignInstance &operator=(const EntityDofRealignInstance &other_permutation);

  /// Destructor
  ~EntityDofRealignInstance();

  /// Return the interpolation point set type that
  /// this reference element represents
  PointSetTag type_id() const;

  /// Return the size of the permutation vector
  Uint size() const;

  /// Return the i-th value
  Uint vertex(const Uint i) const;

  /// Return the code of this permutation
  EntityRealignCode code() const;

  /// Print the internal data
  void print() const;

  private:
  /// Type of reference topology
  PointSetTag m_std_reg_tag;
  EntityRealignCode m_realign_code;
  std::vector<Uint> m_permutation;
};

// ----------------------------------------------------------------------------

inline PointSetTag EntityDofRealignInstance::type_id() const
{
  return m_std_reg_tag;
}

// ----------------------------------------------------------------------------

inline Uint EntityDofRealignInstance::size() const
{
  return m_permutation.size();
}

// ----------------------------------------------------------------------------

inline Uint EntityDofRealignInstance::vertex(const Uint i) const
{
  return m_permutation[i];
}

// ----------------------------------------------------------------------------

inline EntityRealignCode EntityDofRealignInstance::code() const
{
  return m_realign_code;
}

// ----------------------------------------------------------------------------

struct EntityDofRealignFlyweightPolicy
{
  typedef std::pair<PointSetTag, EntityRealignCode> key_type;
};

// ----------------------------------------------------------------------------

typedef common::Flyweight<EntityDofRealignInstance, EntityDofRealignFlyweightPolicy>
    EntityDofRealign;

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
