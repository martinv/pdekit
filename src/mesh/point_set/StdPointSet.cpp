#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/point_set/StdPointSetFactory.hpp"

namespace pdekit
{

namespace mesh
{

namespace detail
{

// ----------------------------------------------------------------------------

const StdPointSetInstance::tag_type StdPointSetInstance::undefined =
    StdPointSetInstance::tag_type(ElemShape::Undefined, P0, PointSetID::Undefined);

// ----------------------------------------------------------------------------

void StdPointSetInstance::construct(const mesh::PointSetTag &tag, StdPointSetInstance &qinstance)
{
  qinstance.m_tag = tag;

  StdPointSetFactory::instance_type &qfactory = StdPointSetFactory::instance();
  const StdPointSetFactory::instance_type::const_product_base_ptr quadrature = qfactory.create(tag);

  qinstance.m_dim   = quadrature->dim();
  qinstance.m_codim = quadrature->codim();
  qinstance.m_coords.resize(quadrature->nb_local_entities());
  qinstance.m_weights.resize(quadrature->nb_local_entities());

  for (Uint e = 0; e < quadrature->nb_local_entities(); ++e)
  {
    quadrature->reference_coords(qinstance.m_coords[e], e);
    quadrature->weights(qinstance.m_weights[e], e);
  }
}

// ----------------------------------------------------------------------------

StdPointSetInstance::StdPointSetInstance()
{
  m_dim   = _0D;
  m_codim = _0D;
  m_coords.resize(0);
  m_weights.resize(0);
}

// ----------------------------------------------------------------------------

StdPointSetInstance::~StdPointSetInstance()
{
}

// ----------------------------------------------------------------------------

const mesh::PointSetTag &StdPointSetInstance::tag() const
{
  return m_tag;
}

// ----------------------------------------------------------------------------

Uint StdPointSetInstance::dim() const
{
  return m_dim;
}

// ----------------------------------------------------------------------------

Uint StdPointSetInstance::codim() const
{
  return m_codim;
}

// ----------------------------------------------------------------------------

Uint StdPointSetInstance::nb_local_entities() const
{
  return m_coords.size();
}

// ----------------------------------------------------------------------------

Uint StdPointSetInstance::size(const Uint local_idx) const
{
  return m_coords[local_idx].rows();
}

// ----------------------------------------------------------------------------

const math::DenseDMat<Real> &StdPointSetInstance::coordinates(const Uint local_idx) const
{
  return m_coords[local_idx];
}

// ----------------------------------------------------------------------------

const math::DenseDVec<Real> &StdPointSetInstance::weights(const Uint local_idx) const
{
  return m_weights[local_idx];
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace mesh

} // namespace pdekit
