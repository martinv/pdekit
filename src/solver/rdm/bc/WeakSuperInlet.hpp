#ifndef PDEKIT_RDM_Weak_Inlet_hpp
#define PDEKIT_RDM_Weak_Inlet_hpp

#include "solver/rdm/bc/WeakBC.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename MeshConfig, typename Physics, Uint BcDim = Physics::DIM - 1>
class WeakSuperInlet
    : public WeakBC<MeshConfig, Physics, BcDim, WeakSuperInlet<MeshConfig, Physics, BcDim>>
{
  private:
  using base = WeakBC<MeshConfig, Physics, BcDim, WeakSuperInlet<MeshConfig, Physics, BcDim>>;
  using geo_cell_metric_type  = typename base::geo_cell_metric_type;
  using sol_cell_metric_type  = typename base::sol_cell_metric_type;
  using flux_cell_metric_type = typename base::flux_cell_metric_type;

  public:
  /// Constructor
  WeakSuperInlet();

  /// Constructor
  WeakSuperInlet(const std::string &name);

  /// Destructor
  ~WeakSuperInlet() override;

  /// Set reference values
  // void set_values(const math::DynamicVector<Real>& state);

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

  /// Values of flux Jacobian, it's inverse and the diagonal
  /// matrix of Jacobian eigenvalues
  typename Physics::JM Rv;
  typename Physics::JM Lv;

  /// Matrix of negative eigenvalues
  typename Physics::JM Dv;
};

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSuperInlet<MeshConfig, Physics, BcDim>::WeakSuperInlet() : base()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSuperInlet<MeshConfig, Physics, BcDim>::WeakSuperInlet(const std::string &name) : base(name)
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
WeakSuperInlet<MeshConfig, Physics, BcDim>::~WeakSuperInlet()
{
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakSuperInlet<MeshConfig, Physics, BcDim>::compute_on_element_weak(
    const Uint idx_in_metric, geo_cell_metric_type const &cell_geo_metric,
    sol_cell_metric_type const &cell_sol_metric, flux_cell_metric_type const &cell_flux_metric)
{
  const Uint nb_qd_pts = cell_geo_metric.nb_qd_pts();
  const Uint nb_nodes  = cell_sol_metric.nb_dof_in_cell();

  // const math::ConstVectorBlock<Real> jacobians = cell_geo_metric.jdet();
  // const math::DynamicVector<Real> &weights = cell_geo_metric.pt_weights();
  const math::DenseConstMatView<Real> qd_pts_phys_coord = cell_geo_metric.interpolated_coords();
  const math::DenseConstMatView<Real> normals_at_qd_pts = cell_geo_metric.normals();

  // const math::DynamicMatrix<Real> &V =
  // cell_geo_metric.reference_sf_values();

  const math::DenseConstMatView<Real> sol_at_qd_pts = cell_sol_metric.field_values();

  // base::m_bface_residual.fill(0.0);

  // Integrate the corrective residual along the boundary face
  // Compute the update coefficients

  math::DenseDVec<Real> &bface_res =
      *(base::m_bface_residual.std_region_data(cell_sol_metric.std_region_type()));

  math::DenseDVec<Real> &bface_update_coeff_vector =
      *(base::m_bface_elem_update_coeff.std_region_data(cell_sol_metric.std_region_type()));
  bface_update_coeff_vector.fill(0.0);

  Real eig_max = 0.0;

  for (Uint q = 0; q < nb_qd_pts; ++q)
  {
    // const Real wq = jacobians[q] * weights[q];

    const math::DenseConstVecView<Real> state_at_qd_pt = sol_at_qd_pts.row_transpose(q);

    Physics::compute_properties(qd_pts_phys_coord.row_transpose(q), state_at_qd_pt,
                                m_solution_gradient, m_phys_properties);

    // std::cout << "Normal = " << base::m_normals[q] << std::endl;

    const math::DenseConstVecView<Real> normal = normals_at_qd_pts.row_transpose(q);

    Physics::flux_jacobian_eigen_structure(m_phys_properties, normal, Rv, Lv, Dv);

    eig_max = 1.e-9;
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      eig_max = std::max(std::abs(Dv(eq, eq)), eig_max);
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

  } // Loop over quadrature points

  //  for (Uint n = 0; n < nb_nodes; ++n)
  //  {
  //    /*
  //    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  //    {
  //      bface_res[n][eq] = 0.0;
  //    }
  //    */
  //    bface_res[n].fill(0.0);
  //  }

  bface_res.fill(0.0);
  base::m_total_facet_residual.fill(0.0);
}

// ----------------------------------------------------------------------------

template <typename MeshConfig, typename Physics, Uint BcDim>
void WeakSuperInlet<MeshConfig, Physics, BcDim>::print_bc_parameters() const
{
  /*
  std::cout << "Supersonic inlet boundary condition:" << std::endl;
  std::cout << "Reference values at infinity = " << m_state_at_infty <<
  std::endl;

  // Square of velocity
  Real vel2 = 0.0;
  for(Uint d = 0; d < Physics::DIM; ++d)
  {
    vel2 += m_state_at_infty[d+1]*m_state_at_infty[d+1];
  }

  const Real p = (1.4-1.)*(m_state_at_infty[Physics::DIM+1] -
  0.5*(vel2)/m_state_at_infty[0]);

  const Real a = std::sqrt(1.4*p/m_state_at_infty[0]);

  const Real Ma = std::sqrt( vel2 )/(m_state_at_infty[0]*a);

  std::cout << "Speed of sound at infinity = " << a << std::endl;
  std::cout << "Mach number at infinity = " << Ma << std::endl;
  */
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
