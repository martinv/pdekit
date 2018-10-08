#include "mesh/EntityDofRealign.hpp"
#include "mesh/std_region/StdRegionFactory.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Class 'EntityRealign'
// ----------------------------------------------------------------------------

const std::pair<PointSetTag, EntityRealignCode> EntityDofRealignInstance::undefined = //
    std::pair<PointSetTag, EntityRealignCode>(
        PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), EntityRealignCode());
// std::pair<Uint,std::string>(StdRegionTypeToID<ElemShapeID::Undefined,P0,PointSetID::Undefined>::value,"");

// ----------------------------------------------------------------------------

void EntityDofRealignInstance::construct(const std::pair<PointSetTag, EntityRealignCode> &key,
                                         EntityDofRealignInstance &entity_permutation)
{
  // Get the interpolation point set from factory to fill the local
  // connectivity and coordinates
  StdRegionFactory::instance_type &factory                         = StdRegionFactory::instance();
  const StdRegionFactory::instance_type::const_product_base_ptr rt = factory.create(key.first);

  entity_permutation.m_std_reg_tag  = key.first;
  entity_permutation.m_realign_code = key.second;

  rt->permutation(key.second, entity_permutation.m_permutation);
}

// ----------------------------------------------------------------------------

EntityDofRealignInstance::EntityDofRealignInstance()
    : m_std_reg_tag(PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined))
{
  m_realign_code.reset();
}

// ----------------------------------------------------------------------------

EntityDofRealignInstance::EntityDofRealignInstance(
    const std::pair<PointSetTag, EntityRealignCode> &key)
{
  construct(key, *this);
}

// ----------------------------------------------------------------------------

EntityDofRealignInstance::EntityDofRealignInstance(
    const EntityDofRealignInstance &other_permutation)
{
  m_std_reg_tag  = other_permutation.m_std_reg_tag;
  m_realign_code = other_permutation.m_realign_code;
  m_permutation.resize(other_permutation.m_permutation.size());
  std::copy(other_permutation.m_permutation.begin(), other_permutation.m_permutation.end(),
            m_permutation.begin());
}

// ----------------------------------------------------------------------------

EntityDofRealignInstance &EntityDofRealignInstance::operator=(
    const EntityDofRealignInstance &other_permutation)
{
  m_std_reg_tag  = other_permutation.m_std_reg_tag;
  m_realign_code = other_permutation.m_realign_code;
  m_permutation.resize(other_permutation.m_permutation.size());
  std::copy(other_permutation.m_permutation.begin(), other_permutation.m_permutation.end(),
            m_permutation.begin());
  return *this;
}

// ----------------------------------------------------------------------------

EntityDofRealignInstance::~EntityDofRealignInstance()
{
}

// ----------------------------------------------------------------------------

void EntityDofRealignInstance::print() const
{
  for (Uint i = 0; i < m_permutation.size(); ++i)
  {
    std::cout << m_permutation[i] << " ";
  }

  std::cout << std::endl;

  std::cout << "Std region of permutation = " << m_std_reg_tag.as_string() << std::endl;
  std::cout << "Permutation code = " << m_realign_code << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
