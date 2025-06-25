#include "mesh/point_set/DofPointSetQuadWarpblend.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

DofPointSetQuadWarpblend::DofPointSetQuadWarpblend(const std::tuple<Uint> &p_order)
    : m_poly_order(std::get<0>(p_order))
{
}

// ----------------------------------------------------------------------------

Uint DofPointSetQuadWarpblend::order() const
{
  return m_poly_order;
}

// ----------------------------------------------------------------------------

Uint DofPointSetQuadWarpblend::dim() const
{
  return std_reg_t::TopoDim;
}

// ----------------------------------------------------------------------------

Uint DofPointSetQuadWarpblend::codim() const
{
  return 0;
}

// ----------------------------------------------------------------------------

Uint DofPointSetQuadWarpblend::nb_local_entities() const
{
  return 1;
}

// ----------------------------------------------------------------------------

Uint DofPointSetQuadWarpblend::size(const Uint local_idx) const
{
  return std_reg_t::nb_dof(m_poly_order);
}

// ----------------------------------------------------------------------------

void DofPointSetQuadWarpblend::reference_coords(math::DenseDMat<Real> &coords,
                                                const Uint local_idx) const
{
  return std_reg_t::fill_coordinates(m_poly_order, coords);
}

// ----------------------------------------------------------------------------

void DofPointSetQuadWarpblend::weights(math::DenseDVec<Real> &wgt, const Uint local_idx) const
{
  wgt.resize(std_reg_t::nb_dof(m_poly_order));
  const Real one_weight = std_reg_t::ref_measure / std_reg_t::nb_dof(m_poly_order);
  for (Uint i = 0; i < wgt.size(); ++i)
  {
    wgt[i] = one_weight;
  }
}

// ----------------------------------------------------------------------------

void DofPointSetQuadWarpblend::permutation(const Uint local_id,
                                           const mesh::EntityRealignCode &permutation_code,
                                           std::vector<Uint> &permutation_vec)
{
  std_reg_t::fill_permutation(m_poly_order, permutation_code, permutation_vec);
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit
