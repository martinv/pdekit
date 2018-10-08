#ifndef PDEKIT_RDM_Weak_Sub_Outlet_BC_hpp
#define PDEKIT_RDM_Weak_Sub_Outlet_BC_hpp

#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakSubOutlet
    : public WeakBC<MeshConfig, Physics, BcDim, WeakSubOutlet<MeshConfig, Physics, BcDim>>
{

  private:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakSubOutlet<MeshConfig, Physics, BcDim>>;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  public:
  /// Default constructor
  WeakSubOutlet();

  /// Constructor
  WeakSubOutlet(const std::string &name);

  /// Destructor
  ~WeakSubOutlet() override;

  /// Set a parameter for the boundary conditions
  void set_parameter(const std::string &param_name, const Real value) override;

  /// Compute the corrective residuals on one element
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  private:
  typename Physics::Properties m_phys_properties;
  typename Physics::Properties::SolGradM m_solution_gradient;
  typename Physics::Properties::FluxV m_ghost_state;
  typename Physics::Properties::FluxV m_ghost_flux_vector;

  // Values of flux Jacobian, it's inverse and the diagonal
  // matrix of Jacobian eigenvalues
  typename Physics::JM Rv;
  typename Physics::JM Lv;
  // Matrix of negative eigenvalues
  typename Physics::JM Dvm;

  // Re-composed Jacobian matrix
  typename Physics::JM GhostJacobian;

  Real m_p_tot;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSubOutlet<MeshConfig, Physics, BcDim>::WeakSubOutlet() : base()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSubOutlet<MeshConfig, Physics, BcDim>::WeakSubOutlet(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSubOutlet<MeshConfig, Physics, BcDim>::~WeakSubOutlet()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakSubOutlet<MeshConfig, Physics, BcDim>::set_parameter(const std::string &param_name,
                                                              const Real value)
{
  if (param_name == "total_pressure")
  {
    m_p_tot = value;
  }
  else
  {
    std::cerr << "WeakSubOutlet::set_parameter: unknown parameter name '" << param_name << "'"
              << std::endl;
  }
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakSubOutlet<MeshConfig, Physics, BcDim>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
  const Uint nb_qd_pts = cell_geo_metric.nb_qd_pts();
  const Uint nb_nodes  = cell_sol_metric.nb_dof_in_cell();

  const math::DenseConstVecView<Real> jacobians         = cell_geo_metric.jdet();
  const math::DenseDVec<Real> &weights                  = cell_geo_metric.pt_weights();
  const math::DenseConstMatView<Real> qd_pts_phys_coord = cell_geo_metric.interpolated_coords();
  const math::DenseConstMatView<Real> normals_at_qd_pts = cell_geo_metric.normals();

  const math::DenseDMat<Real> &V = cell_sol_metric.reference_sf_values();

  const math::DenseConstMatView<Real> sol_at_qd_pts = cell_sol_metric.field_values();

  math::DenseDVec<Real> &bface_res =
      *(base::m_bface_residual.std_region_data(cell_sol_metric.std_region_type()));

  math::DenseDVec<Real> &bface_update_coeff_vector =
      *(base::m_bface_elem_update_coeff.std_region_data(cell_sol_metric.std_region_type()));

  // Integrate the corrective residual along the boundary face
  Real eig_max = 0.0;

  base::m_total_facet_residual.fill(0.0);

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobians[q] * weights[q];

    const math::DenseConstVecView<Real> state_at_qd_pt = sol_at_qd_pts.row_transpose(q);

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), state_at_qd_pt,
                                m_solution_gradient, m_phys_properties);

    const math::DenseConstVecView<Real> normal = normals_at_qd_pts.row_transpose(q);

    Physics::flux_jacobian_eigen_structure(m_phys_properties, normal, Rv, Lv, Dvm);

    // Density and momentum terms are just copied into the ghost state:
    for (Uint i = 0; i < (Physics::NEQ - 1); ++i)
    {
      m_ghost_state[i] = state_at_qd_pt[i];
    }

    // This variable is rho^2 *(u^2+v^2+w^2)
    Real rho2u2 = 0.0;
    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      rho2u2 += state_at_qd_pt[d + 1] * state_at_qd_pt[d + 1];
    }

    m_ghost_state[Physics::NEQ - 1] =
        m_p_tot / (m_phys_properties.gamma_minus_1) + 0.5 * (rho2u2) / state_at_qd_pt[0];

    // m_ghost_state(0,3) = base::m_uq(q,3);

    eig_max = 1.e-9;

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      eig_max     = std::max(std::abs(Dvm(eq, eq)), eig_max);
      Dvm(eq, eq) = std::min(Dvm(eq, eq), 0.0);
      m_ghost_state[eq] -= state_at_qd_pt[eq];
    }

    /*
    for (Uint n = 0; n < bface_res.size(); ++n)
    {
      bface_update_coeff_vector[n] += wq * eig_max;
    }
    */

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      bface_update_coeff_vector[n] = std::max(bface_update_coeff_vector[n], eig_max);
    }

    GhostJacobian.fill(0.0);
    GhostJacobian = Rv * Dvm * Lv;

    m_ghost_flux_vector = GhostJacobian * m_ghost_state;
    base::m_total_facet_residual += wq * m_ghost_flux_vector;

    for (Uint n = 0; n < nb_nodes; ++n)
    {
      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        bface_res[n * Physics::NEQ + eq] += V(q, n) * wq * m_ghost_flux_vector[eq];
      }
    }
  } // Loop over integration points
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
