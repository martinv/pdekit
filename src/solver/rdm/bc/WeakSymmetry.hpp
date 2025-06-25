#ifndef PDEKIT_RDM_Weak_Symmetry_hpp
#define PDEKIT_RDM_Weak_Symmetry_hpp

#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakSymmetry
    : public WeakBC<MeshConfig, Physics, BcDim, WeakSymmetry<MeshConfig, Physics, BcDim>>
{
  private:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakSymmetry<MeshConfig, Physics, BcDim>>;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  public:
  /// Default constructor
  WeakSymmetry();

  /// Constructor
  WeakSymmetry(const std::string &name);

  /// Destructor
  ~WeakSymmetry() override;

  /// Compute the corrective residuals on one element
  void compute_on_element_weak(const Uint idx_in_metric,
                               geo_cell_metric_type const &cell_geo_metric,
                               sol_cell_metric_type const &cell_sol_metric,
                               flux_cell_metric_type const &cell_flux_metric);

  private:
  typename Physics::Properties m_phys_props;
  typename Physics::Properties::SolGradM m_solution_gradient;
  typename Physics::Properties::FluxV m_ghost_state;
  typename Physics::Properties::FluxV m_ghost_flux_vector;

  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename Physics::JM Rv;
  typename Physics::JM Lv;
  /// Matrix of negative eigenvalues
  typename Physics::JM Dvm;

  /// Re-composed Jacobian matrix
  typename Physics::JM GhostJacobian;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSymmetry<MeshConfig, Physics, BcDim>::WeakSymmetry() : base()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSymmetry<MeshConfig, Physics, BcDim>::WeakSymmetry(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSymmetry<MeshConfig, Physics, BcDim>::~WeakSymmetry()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakSymmetry<MeshConfig, Physics, BcDim>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
  // base::m_bface_residual.fill(0.0);

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

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), sol_at_qd_pts.row_transpose(q),
                                m_solution_gradient, m_phys_props);

    const math::DenseConstVecView<Real> normal = normals_at_qd_pts.row_transpose(q);

    Physics::flux_jacobian_eigen_structure(m_phys_props, normal, Rv, Lv, Dvm);

    // Reflected velocity (flipped with respect to boundary normal)
    // If original velocity vector was v and the unit boundary normal is n,
    // then the flipped velocity will be v_flip = v - 2*(v,n) n

    // Scalar product v . n :
    Real v_dot_n = 0.0;

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      v_dot_n += normal[d] * m_phys_props.V[d];
    }

    m_ghost_state[0]                = state_at_qd_pt[0];
    m_ghost_state[Physics::NEQ - 1] = state_at_qd_pt[Physics::NEQ - 1];

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      // Momentum terms:      density *        (              reflected
      // velocity
      // )
      m_ghost_state[d + 1] = state_at_qd_pt[0] * (m_phys_props.V[d] - 2. * v_dot_n * normal[d]);
    }

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
