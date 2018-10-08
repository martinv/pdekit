#include "mesh/point_set/QuadraturePermutation.hpp"
#include "mesh/point_set/StdPointSetFactory.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------
// Class 'QuadraturePermutationInstance' - static variables
// ----------------------------------------------------------------------------

const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode>
    QuadraturePermutationInstance::undefined =
        std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode>(
            mesh::PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined), 0,
            mesh::EntityRealignCode());

// ----------------------------------------------------------------------------

void QuadraturePermutationInstance::construct(
    const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> &key,
    QuadraturePermutationInstance &quad_permutation)
{
  // Get the interpolation point set from factory to fill the local
  // connectivity and coordinates
  StdPointSetFactory::instance_type &factory = StdPointSetFactory::instance();
  const StdPointSetFactory::instance_type::const_product_base_ptr quad =
      factory.create(std::get<0>(key));

  quad_permutation.m_quad_tag         = std::get<0>(key);
  quad_permutation.m_local_id         = std::get<1>(key);
  quad_permutation.m_permutation_code = std::get<2>(key);

  quad->permutation(std::get<1>(key), std::get<2>(key), quad_permutation.m_permutation);
}

// ----------------------------------------------------------------------------

QuadraturePermutationInstance::QuadraturePermutationInstance()
    : m_quad_tag(mesh::PointSetTag(ElemShape::Undefined, P0, PointSetID::Undefined)), m_local_id(0)
{
  m_permutation_code.reset();
}

// ----------------------------------------------------------------------------

QuadraturePermutationInstance::QuadraturePermutationInstance(
    const std::tuple<mesh::PointSetTag, Uint, mesh::EntityRealignCode> &key)
{
  construct(key, *this);
}

// ----------------------------------------------------------------------------

QuadraturePermutationInstance::QuadraturePermutationInstance(
    const QuadraturePermutationInstance &other_permutation)
{
  m_quad_tag         = other_permutation.m_quad_tag;
  m_local_id         = other_permutation.m_local_id;
  m_permutation_code = other_permutation.m_permutation_code;
  m_permutation.resize(other_permutation.m_permutation.size());
  std::copy(other_permutation.m_permutation.begin(), other_permutation.m_permutation.end(),
            m_permutation.begin());
}

// ----------------------------------------------------------------------------

QuadraturePermutationInstance &QuadraturePermutationInstance::operator=(
    const QuadraturePermutationInstance &other_permutation)
{
  m_quad_tag         = other_permutation.m_quad_tag;
  m_permutation_code = other_permutation.m_permutation_code;
  m_permutation.resize(other_permutation.m_permutation.size());
  std::copy(other_permutation.m_permutation.begin(), other_permutation.m_permutation.end(),
            m_permutation.begin());
  return *this;
}

// ----------------------------------------------------------------------------

QuadraturePermutationInstance::~QuadraturePermutationInstance()
{
}

// ----------------------------------------------------------------------------

void QuadraturePermutationInstance::print() const
{
  for (Uint i = 0; i < m_permutation.size(); ++i)
  {
    std::cout << m_permutation[i] << " ";
  }

  std::cout << std::endl;

  std::cout << "Std region of permutation = " << m_quad_tag.as_string() << std::endl;
  std::cout << "Permutation code = " << m_permutation_code << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
