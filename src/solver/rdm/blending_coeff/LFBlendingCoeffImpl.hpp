#ifndef PDEKIT_Solver_RDM_Lax_Friedrichs_Blending_Coeff_hpp
#define PDEKIT_Solver_RDM_Lax_Friedrichs_Blending_Coeff_hpp

#include "math/DenseDVec.hpp"
#include "math/MathConstants.hpp"
#include "mesh/MeshConfig.hpp"
#include "solver/rdm/blending_coeff/BlendingCoeffImplBase.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

namespace detail
{

template <typename MeshConfig>
class LFBlendingCoeffImpl : public BlendingCoeffImplBase<MeshConfig>
{
  public:
  /// Default constructor
  LFBlendingCoeffImpl();

  /// Default destructor
  ~LFBlendingCoeffImpl() override;

  /// Resize internal data
  void setup(typename result_of::dof_map_t<MeshConfig> const &sol_dofs) override;

  /// Evaluate the blending coefficient
  void calculate_impl(typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
                      const interpolation::VectorMeshFunction<Real> &u,
                      interpolation::ScalarMeshFunction<Real> &function) override;

  /// Clear all internal data
  void clear() override;

  /// Set a parameter
  void set_param(const std::string &param_name, const Real param_value) override;

  private:
  /// DATA
  std::vector<Real> m_theta_nodal;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
LFBlendingCoeffImpl<MeshConfig>::LFBlendingCoeffImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
LFBlendingCoeffImpl<MeshConfig>::~LFBlendingCoeffImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LFBlendingCoeffImpl<MeshConfig>::setup(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs)
{
  m_theta_nodal.resize(sol_dofs.nb_nodes());
  m_theta_nodal.assign(sol_dofs.nb_nodes(), 1.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LFBlendingCoeffImpl<MeshConfig>::calculate_impl(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &function)
{

  typedef interpolation::VectorMeshFunction<Real>::const_entry_type node_value_type;

  m_theta_nodal.assign(m_theta_nodal.size(), 0.0);

  const Real eps = 1.e-12;
  Real u_avg, ksi_elem;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_dofs.active_cell(mesh::ActiveIdx(c));

    // Compute average state u_avg in element
    u_avg = 0.0;
    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const node_value_type u_node = u.const_value(cell.vertex(v));
      u_avg += u_node[0];
    }
    u_avg /= cell.nb_vert();

    // Find the value max (|u-u_avg|/(|u|+|u_avg| + eps) for all u in one
    // element
    ksi_elem = 0.0;
    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const node_value_type u_node = u.const_value(cell.vertex(v));
      ksi_elem                     = std::max(ksi_elem, std::abs(u_node[0] - u_avg) /
                                        (std::abs(u_node[0]) + std::abs(u_avg) + eps));
    }

    // Update all theta values in element so that they are at least
    // as big as currently computed 'ksi_elem'
    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      m_theta_nodal[cell.vertex(v)] = std::max(m_theta_nodal[cell.vertex(v)], ksi_elem);
    }
  } // Loop over cells

  // m_theta = m_theta_nodal;

  Real theta_elem;

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_dofs.active_cell(mesh::ActiveIdx(c));

    theta_elem = 0.0;
    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      theta_elem = std::max(theta_elem, m_theta_nodal[cell.vertex(v)]);
    }

    function[c] = 1. - theta_elem;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LFBlendingCoeffImpl<MeshConfig>::clear()
{
  m_theta_nodal.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void LFBlendingCoeffImpl<MeshConfig>::set_param(const std::string &param_name,
                                                const Real param_value)
{
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
