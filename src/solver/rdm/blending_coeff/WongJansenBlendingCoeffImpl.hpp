#ifndef PDEKIT_Solver_RDM_Wong_Jansen_Blending_Coeff_hpp
#define PDEKIT_Solver_RDM_Wong_Jansen_Blending_Coeff_hpp

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
class WongJansenBlendingCoeffImpl : public BlendingCoeffImplBase<MeshConfig>
{
  public:
  /// Default constructor
  WongJansenBlendingCoeffImpl();

  /// Default destructor
  ~WongJansenBlendingCoeffImpl() override;

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

  /// Freestream solution/velocity value
  Real m_u_inf;

  std::vector<Real> m_theta_nodal;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig>
WongJansenBlendingCoeffImpl<MeshConfig>::WongJansenBlendingCoeffImpl() : m_u_inf(1.e-14)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
WongJansenBlendingCoeffImpl<MeshConfig>::~WongJansenBlendingCoeffImpl()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void WongJansenBlendingCoeffImpl<MeshConfig>::setup(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs)
{
  m_theta_nodal.resize(sol_dofs.nb_nodes());
  m_theta_nodal.assign(sol_dofs.nb_nodes(), 1.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void WongJansenBlendingCoeffImpl<MeshConfig>::calculate_impl(
    typename result_of::dof_map_t<MeshConfig> const &sol_dofs,
    const interpolation::VectorMeshFunction<Real> &u,
    interpolation::ScalarMeshFunction<Real> &function)
{
  // const Real eps = 1.e-1;
  const Real eps    = 1.0;
  const Real u_inf4 = std::pow(m_u_inf, 4);

  for (Uint c = 0; c < sol_dofs.nb_active_cells(); ++c)
  {
    const mesh::MeshEntity cell = sol_dofs.active_cell(mesh::ActiveIdx(c));

    Real nom   = 0.0;
    Real denom = eps * u_inf4;
    Real u_avg = 0.0;

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const interpolation::VectorMeshFunction<Real>::const_entry_type node_value =
          u.const_value(cell.vertex(v));
      u_avg += node_value[0];
    }

    u_avg /= cell.nb_vert();

    for (Uint v = 0; v < cell.nb_vert(); ++v)
    {
      const interpolation::VectorMeshFunction<Real>::const_entry_type node_value =
          u.const_value(cell.vertex(v));
      const Real u_i = node_value[0];

      const Real term = (u_i * u_i - u_avg * u_avg) * (u_i * u_i - u_avg * u_avg);
      nom += term;
      denom += term;
    }

    function[cell.idx()] = nom / denom;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void WongJansenBlendingCoeffImpl<MeshConfig>::clear()
{
  m_theta_nodal.resize(0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig>
void WongJansenBlendingCoeffImpl<MeshConfig>::set_param(const std::string &param_name,
                                                        const Real param_value)
{
  if (param_name == "u_inf")
  {
    m_u_inf = param_value;
  }
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
