#ifndef PDEKIT_RDM_Weak_Farfield_hpp
#define PDEKIT_RDM_Weak_Farfield_hpp

#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakFarfield
    : public WeakBC<MeshConfig, Physics, BcDim, WeakFarfield<MeshConfig, Physics, BcDim>>
{

  private:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakFarfield<MeshConfig, Physics, BcDim>>;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  public:
  /// Constructor
  WeakFarfield();

  /// Constructor
  WeakFarfield(const std::string &name);

  /// Destructor
  ~WeakFarfield() override;

  /// Set reference values
  void set_reference_state(const math::DenseDVec<Real> &state) override;

  /// Compute the corrective residuals on one element
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  /// Print the reference values stored in this bc
  void print_bc_parameters() const override;

  private:
  typename Physics::Properties m_phys_properties;
  typename Physics::Properties::SolGradM m_solution_gradient;
  typename Physics::Properties::FluxV m_ghost_state;
  typename Physics::Properties::FluxV m_reference_state;
  typename Physics::Properties::FluxV m_ghost_flux_vector;

  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename Physics::JM Rv;
  typename Physics::JM Lv;

  /// Matrix of negative eigenvalues
  typename Physics::JM Dvm;

  /// Re-composed Jacobian matrix
  typename Physics::JM GhostJacobian;

  /// Pressure corresponding to reference state
  Real m_p_ref;

  /// Velocity square corresponding to reference state
  Real m_v2_ref;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakFarfield<MeshConfig, Physics, BcDim>::WeakFarfield() : base()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakFarfield<MeshConfig, Physics, BcDim>::WeakFarfield(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakFarfield<MeshConfig, Physics, BcDim>::~WeakFarfield()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakFarfield<MeshConfig, Physics, BcDim>::set_reference_state(
    const math::DenseDVec<Real> &state)
{
  for (Uint i = 0; i < Physics::NEQ; ++i)
  {
    m_reference_state[i] = state[i];
  }

  Real v2 = 0;
  for (Uint d = 1; d < (Physics::NEQ - 1); ++d)
  {
    v2 += state[d] * state[d];
  }
  v2 /= (state[0] * state[0]);

  m_v2_ref = v2;

  m_p_ref = (1.4 - 1.) * (state[Physics::NEQ - 1] - 0.5 * state[0] * v2);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakFarfield<MeshConfig, Physics, BcDim>::compute_on_element_weak(
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

  // base::m_bface_residual.fill(0.0);
  // bface_update_coeff_vector.fill(0.0);

  // Integrate the corrective residual along the boundary face
  Real eig_max = 1.e-8;

  base::m_total_facet_residual.fill(0.0);

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    const Real wq = jacobians[q] * weights[q];

    const math::DenseConstVecView<Real> state_at_qd_pt = sol_at_qd_pts.row_transpose(q);

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), state_at_qd_pt,
                                m_solution_gradient, m_phys_properties);

    // std::cout << "Normal = " << base::m_normals[q] << std::endl;

    const math::DenseConstVecView<Real> normal = normals_at_qd_pts.row_transpose(q);

    Physics::flux_jacobian_eigen_structure(m_phys_properties, normal, Rv, Lv, Dvm);

    // const Real e_ref = m_p_ref/(1.4-1.0) + 0.5 * state_at_qd_pt[0] *
    // m_phys_properties.v2;

    // Reference state is the same as the state at quadrature point, but has
    // imposed pressure
    // We set energy which is computed from the reference pressure
    // m_state_at_infty = state_at_qd_pt;
    // m_state_at_infty[Physics::NEQ-1] = e_ref;

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      eig_max     = std::max(std::abs(Dvm(eq, eq)), eig_max);
      Dvm(eq, eq) = std::min(Dvm(eq, eq), 0.0);
      // m_ghost_state[eq] = m_reference_state[eq] - state_at_qd_pt[eq];
    }

    m_ghost_state = m_reference_state;

#if 0
    Real vn = 0.0; // Normal velocity
    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      vn += m_phys_properties.V[d] * normal[d];
    }

    if (vn < 0.0) // Inlet
    {
      // Prescribe the first n-1 variables (all except energy/pressure)
      for (Uint eq = 0; eq < (Physics::NEQ - 1); ++eq)
      {
        m_ghost_state[eq] = m_reference_state[eq];
      }

      // Subsonic inlet - extrapolate pressure
      if (std::abs(vn) < m_phys_properties.a)
      {
        const Real rhoE = m_phys_properties.P / (1.4 - 1.) + 0.5 * m_reference_state[0] * m_v2_ref;
        m_ghost_state[Physics::NEQ - 1] = rhoE;
      }
      // Supersonic inlet - presribe everything
      else
      {
        m_ghost_state[Physics::NEQ - 1] = m_reference_state[Physics::NEQ - 1];
      }
    }
    // Outlet
    else
    {
      // Extrapolate the first n-1 variables (all except energy/pressure)
      for (Uint eq = 0; eq < (Physics::NEQ - 1); ++eq)
      {
        m_ghost_state[eq] = state_at_qd_pt[eq];
      }

      // Subsonic outlet - prescribe pressure
      if (std::abs(vn) < m_phys_properties.a)
      {
        const Real rhoE = m_p_ref / (1.4 - 1.) + 0.5 * m_phys_properties.rho * m_phys_properties.v2;
        m_ghost_state[Physics::NEQ - 1] = rhoE;
      }
      // Supersonic outlet - extrapolated everything
      else
      {
        m_ghost_state[Physics::NEQ - 1] = state_at_qd_pt[Physics::NEQ - 1];
      }
    }
#else

    // Now use m_ghost_state to store the value of (u_in/out - u_q )

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      m_ghost_state[eq] -= state_at_qd_pt[eq];
    }
#endif

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

    /*
    std::cout << "BFACE_RES: " << bface_res.size() << std::endl;
    std::cout << "n = [0-" << nb_nodes-1 << "]" << std::endl;
    std::cout << "V: " << V.rows() << " x " << V.cols() << std::endl;
    std::cout << "Nb. qd. pts = " << nb_qd_pts << std::endl;
    std::cout << "Ghost flux vector = " << m_ghost_flux_vector.size() <<
    std::endl;
    */

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

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakFarfield<MeshConfig, Physics, BcDim>::print_bc_parameters() const
{
  std::cout << "WeakFarfield:: reference pressure = " << m_p_ref << std::endl;
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
