#ifndef PDEKIT_Solver_RDM_Lax_Friedrichs_hpp
#define PDEKIT_Solver_RDM_Lax_Friedrichs_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"
#include "solver/rdm/cellsplitters/LFResidualLimiter.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Tag class to mark Lax-Friedrichs scheme

// ----------------------------------------------------------------------------

class LF
{
};

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics>
class LFImplementation
{
  public:
  enum
  {
    is_var_beta_type_rdm = 0
  };

  /// TYPEDEFS

  typedef Physics phys_model;
  typedef PGRDLFMethodData<Physics> method_data;

  /// Constructor
  LFImplementation();

  /// Destructor
  ~LFImplementation();

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, typename MD3>
  void compute_adv_residuals(math::DenseConstMatView<Real> const &sol_nodal_values,
                             interpolation::CellGeoMetric<MD1> const &CGM,
                             interpolation::CellSolutionMetric<MD2> const &CSM,
                             interpolation::CellFluxMetric<MD3, Physics> const &CFM,
                             PGRDLFMethodData<Physics> &MD);

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, typename MD3>
  void compute_adv_reaction_residuals(math::DenseConstMatView<Real> const &sol_nodal_values,
                                      interpolation::CellGeoMetric<MD1> const &CGM,
                                      interpolation::CellSolutionMetric<MD2> const &CSM,
                                      interpolation::CellFluxMetric<MD3, Physics> const &CFM,
                                      interpolation::CellSolutionMetric<MD2> const &CSrcM,
                                      PGRDLFMethodData<Physics> &MD);

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2>
  void compute_adv_residuals_no_flux_metric(math::DenseConstMatView<Real> const &sol_nodal_values,
                                            interpolation::CellGeoMetric<MD1> const &CGM,
                                            interpolation::CellSolutionMetric<MD2> const &CSM,
                                            typename Physics::FluxV const &cell_residual,
                                            PGRDLFMethodData<Physics> &MD);

  private:
  /// TYPES
  typedef LFResidualLimiter<Physics, Physics::NEQ> residual_limiter_type;
};

// ----------------------------------------------------------------------------

template <typename Physics>
LFImplementation<Physics>::LFImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
LFImplementation<Physics>::~LFImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3>
void LFImplementation<Physics>::compute_adv_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    interpolation::CellGeoMetric<MD1> const &CGM, interpolation::CellSolutionMetric<MD2> const &CSM,
    interpolation::CellFluxMetric<MD3, Physics> const &CFM, PGRDLFMethodData<Physics> &MD)
{
  // Reset the cell residuals
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n].fill(0.0);
    MD.m_LF_proj_vec[n].fill(0.0);
  }

  // Reset the cell update coefficient
  MD.m_k_LF.fill(0.0);
  MD.m_elem_wave_speed.fill(0.0);

  MD.m_u_avg.fill(0.0);

  // sol_nodal_values is a matrix of dimension
  // [(nb. nodes) x (nb. equations)]

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_u_avg += sol_nodal_values.row_transpose(n);
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    MD.m_u_avg[eq] *= 1. / MD.m_nb_nodes;
  }

  // Quadrature weights
  const math::DenseDVec<Real> &w = CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;
  // Solution in quadrature points

  const math::DenseConstMatView<Real> uq = CSM.field_values();

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = CSM.field_derivatives(d);
    MD.m_grad_F[d] = CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // ==============================================================
  // PART I: COMPUTE THE (NONLIMITED) LF RESIDUALS (CENTRAL METHOD)
  // ==============================================================

  MD.m_total_cell_res.fill(0.0);
  MD.m_cell_volume = 0.0;

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);
    // Compute the gradient of each shape function in physical space
    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    const Real wj_q = w[q] * jdet[q];
    MD.m_cell_volume += wj_q;

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      // Accumulate to the projection vectors ( int{ \nabla \varphi_i } )
      MD.m_LF_proj_vec[n] += (wj_q / Physics::DIM) * MD.m_grad_sf_at_pt_phys;

      /*
      Physics::flux_jacobian_eigen_structure(MD.m_props,
      MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv, MD.Dvp);
      */

      Physics::flux_jacobian_eigen_values(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.m_eigenvalues);

      MD.m_max_eigenvalue = 1.e-9;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.m_eigenvalues[e]));
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));

        // THIS DOESN'T WORK (STABILITY PROBLEMS)!
        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      MD.m_k_LF[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      /*
      // For the moment, we use MD.m_elem_update_coeff[n] to store the
      value of
      // integral { \nabla(F) \cdot \nabla \varphi_n  }
      MD.m_elem_update_coeff[n] += wj_q * dot_prod;
      */

    } // Loop over nodes

    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    MD.m_total_cell_res += wj_q * MD.m_res_at_point;

  } // Loop over quadrature points

  // ----------------------------------------
  // Compute LF alpha (stability coefficient)
  // ----------------------------------------

  // At this point, MD.m_elem_wave_speed[n] holds the maximum spectral radius
  // of ( A \cdot \nabla \varphi_n ) sampled at quadrature points

  MD.m_LF_alpha = MD.m_k_LF[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_LF_alpha = std::max(std::abs(MD.m_k_LF[n]), MD.m_LF_alpha);
  }

  // -----------------------------------------
  // Finish computation of update coefficients
  // -----------------------------------------

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= MD.m_cell_volume;
  }

  // ---------------------------
  // Compute nodal residuals
  // ---------------------------
  MD.m_elem_node_res[0] = 1.0 / MD.m_nb_nodes * MD.m_total_cell_res;

  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] = MD.m_elem_node_res[0];
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] += MD.m_LF_alpha * (sol_nodal_values.row_transpose(n) - MD.m_u_avg);
  }

  // ===============================================================
  // PART II: APPLY LIMITING
  // ===============================================================
  residual_limiter_type::limit_cell_residuals(MD.m_total_cell_res, MD.m_sum_x_p, MD);

#if 1
  // ===============================================================
  // PART III: ADD FILTER (SUPG-LIKE) TERM
  // ===============================================================

  // Compute (approximately) the barycenter of element by computing
  // the mean value of coordinate function
  // At the same time compute 'mean' gradient at barycenter
  MD.m_barycenter.fill(0.0);
  MD.m_grad_u_at_point.fill(0.0);

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    const Real wj_q = w[q] * jdet[q];

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      MD.m_barycenter[d] += wj_q * Xq(q, d);

      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        MD.m_grad_u_at_point(eq, d) += wj_q * MD.m_grad_u[d](q, eq);
      }
    }
  }

  const Real inv_cell_volume = 1. / MD.m_cell_volume;

  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    MD.m_barycenter[d] *= inv_cell_volume;
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      MD.m_grad_u_at_point(eq, d) *= inv_cell_volume;
    }
  }

  // Use sum_Kp and inv_sum_Kp to compute and hold
  // the stabilization matrix
  MD.m_grad_u_at_point.fill(0.0);
  Physics::compute_properties(MD.m_barycenter, MD.m_u_avg, MD.m_grad_u_at_point, MD.m_props);

  MD.sum_Kp.fill(0.0);

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_LF_proj_vec[n], MD.Rv, MD.Lv, MD.Dvp);

    for (Uint e = 0; e < Physics::NEQ; ++e)
    {
      MD.Dvp(e, e) = std::max(1.e-9, MD.Dvp(e, e));
    }

    MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;
    MD.sum_Kp += MD.m_Kp[n];
  }

  // MD.inv_sum_Kp will hold the stabilization matrix 'Sigma' in LF filtering
  // term
  MD.sum_Kp.inv(MD.inv_sum_Kp);

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    // Compute the Jacobian matrices of fluxes (stored in MD.m_dFdu)
    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);
    Physics::residual(MD.m_props, MD.m_dFdu, MD.m_res_at_point);

    // Equation residual: div(F) evaluated at quadrature point q
    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Compute the gradient of each shape function in physical space
    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      // Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      MD.m_J_proj.fill(0.0);
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        MD.m_J_proj += MD.m_grad_sf_at_pt_phys[dim] * MD.m_dFdu[dim];
      }

      const Real wj_q = w[q] * jdet[q];

      MD.m_flux_integral_in_elem = MD.m_J_proj * (MD.inv_sum_Kp * MD.m_res_at_point);
      MD.m_elem_node_res[n] +=
          (MD.m_LF_upwind_stab * 0.5 * MD.m_cell_volume * wj_q) * MD.m_flux_integral_in_elem;

    } // Loop over element nodes

  } // Loop over quadrature points
#endif
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3>
void LFImplementation<Physics>::compute_adv_reaction_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    interpolation::CellGeoMetric<MD1> const &CGM, interpolation::CellSolutionMetric<MD2> const &CSM,
    interpolation::CellFluxMetric<MD3, Physics> const &CFM,
    interpolation::CellSolutionMetric<MD2> const &CSrcM, PGRDLFMethodData<Physics> &MD)
{
  // Reset the cell residuals
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n].fill(0.0);
    MD.m_LF_proj_vec[n].fill(0.0);
  }

  // Reset the cell update coefficient
  MD.m_elem_wave_speed.fill(0.0);

  MD.m_u_avg.fill(0.0);

  // sol_nodal_values is a matrix of dimension
  // [(nb. nodes) x (nb. equations)]

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_u_avg += sol_nodal_values.row_transpose(n);
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    MD.m_u_avg[eq] *= 1. / MD.m_nb_nodes;
  }

  // Quadrature weights
  const math::DenseDVec<Real> &w = CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;
  // Solution in quadrature points

  const math::DenseConstMatView<Real> uq = CSM.field_values();

  // Source field in quadrature points
  const math::DenseConstMatView<Real> sq = CSrcM.field_values();

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = CSM.field_derivatives(d);
    MD.m_grad_F[d] = CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // ==============================================================
  // PART I: COMPUTE THE (NONLIMITED) LF RESIDUALS (CENTRAL METHOD)
  // ==============================================================

  MD.m_total_cell_res.fill(0.0);
  MD.m_cell_volume = 0.0;

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);
    // Compute the gradient of each shape function in physical space
    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    const Real wj_q = w[q] * jdet[q];
    MD.m_cell_volume += wj_q;

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      // Accumulate to the projection vectors ( int{ \nabla \varphi_i } )
      MD.m_LF_proj_vec[n] += (wj_q / Physics::DIM) * MD.m_grad_sf_at_pt_phys;

      /*
      Physics::flux_jacobian_eigen_structure(MD.m_props,
      MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv, MD.Dvp);
      */

      Physics::flux_jacobian_eigen_values(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.m_eigenvalues);

      MD.m_max_eigenvalue = 1.e-9;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.m_eigenvalues[e]));
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));

        // THIS DOESN'T WORK (STABILITY PROBLEMS)!
        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      /*
      // For the moment, we use MD.m_elem_update_coeff[n] to store the
      value of
      // integral { \nabla(F) \cdot \nabla \varphi_n  }
      MD.m_elem_update_coeff[n] += wj_q * dot_prod;
      */

    } // Loop over nodes

    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    MD.m_res_at_point -= sq.row_transpose(q);

    MD.m_total_cell_res += wj_q * MD.m_res_at_point;

  } // Loop over quadrature points

  // ----------------------------------------
  // Compute LF alpha (stability coefficient)
  // ----------------------------------------

  // At this point, MD.m_elem_wave_speed[n] holds the maximum spectral radius
  // of ( A \cdot \nabla \varphi_n ) sampled at quadrature points

  MD.m_LF_alpha = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_LF_alpha = std::max(std::abs(MD.m_elem_wave_speed[n]), MD.m_LF_alpha);
  }

  // -----------------------------------------
  // Finish computation of update coefficients
  // -----------------------------------------

  /*
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= MD.m_cell_volume;
  }
  */

  // ---------------------------
  // Compute nodal residuals
  // ---------------------------
  MD.m_elem_node_res[0] = 1.0 / MD.m_nb_nodes * MD.m_total_cell_res;

  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] = MD.m_elem_node_res[0];
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] += MD.m_LF_alpha * (sol_nodal_values.row_transpose(n) - MD.m_u_avg);
  }

  // ===============================================================
  // PART II: APPLY LIMITING
  // ===============================================================
  residual_limiter_type::limit_cell_residuals(MD.m_total_cell_res, MD.m_sum_x_p, MD);

#if 1
  // ===============================================================
  // PART III: ADD FILTER (SUPG-LIKE) TERM
  // ===============================================================

  // Compute (approximately) the barycenter of element by computing
  // the mean value of coordinate function
  // At the same time compute 'mean' gradient at barycenter
  MD.m_barycenter.fill(0.0);
  MD.m_grad_u_at_point.fill(0.0);

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    const Real wj_q = w[q] * jdet[q];

    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      MD.m_barycenter[d] += wj_q * Xq(q, d);

      for (Uint eq = 0; eq < Physics::NEQ; ++eq)
      {
        MD.m_grad_u_at_point(eq, d) += wj_q * MD.m_grad_u[d](q, eq);
      }
    }
  }

  const Real inv_cell_volume = 1. / MD.m_cell_volume;

  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    MD.m_barycenter[d] *= inv_cell_volume;
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      MD.m_grad_u_at_point(eq, d) *= inv_cell_volume;
    }
  }

  // Use sum_Kp and inv_sum_Kp to compute and hold
  // the stabilization matrix
  MD.m_grad_u_at_point.fill(0.0);
  Physics::compute_properties(MD.m_barycenter, MD.m_u_avg, MD.m_grad_u_at_point, MD.m_props);

  MD.sum_Kp.fill(0.0);

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_LF_proj_vec[n], MD.Rv, MD.Lv, MD.Dvp);

    for (Uint e = 0; e < Physics::NEQ; ++e)
    {
      MD.Dvp(e, e) = std::max(1.e-9, MD.Dvp(e, e));
    }

    MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;
    MD.sum_Kp += MD.m_Kp[n];
  }

  // MD.inv_sum_Kp will hold the stabilization matrix 'Sigma' in LF filtering
  // term
  MD.sum_Kp.inv(MD.inv_sum_Kp);

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    // Compute the Jacobian matrices of fluxes (stored in MD.m_dFdu)
    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);
    Physics::residual(MD.m_props, MD.m_dFdu, MD.m_res_at_point);

    // Equation residual: div(F) evaluated at quadrature point q
    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Compute the gradient of each shape function in physical space
    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      // Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      MD.m_J_proj.fill(0.0);
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        MD.m_J_proj += MD.m_grad_sf_at_pt_phys[dim] * MD.m_dFdu[dim];
      }

      const Real wj_q = w[q] * jdet[q];

      MD.m_flux_integral_in_elem = MD.m_J_proj * (MD.inv_sum_Kp * MD.m_res_at_point);
      MD.m_elem_node_res[n] +=
          (MD.m_LF_upwind_stab * 0.5 * MD.m_cell_volume * wj_q) * MD.m_flux_integral_in_elem;

    } // Loop over element nodes

  } // Loop over quadrature points
#endif
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2>
void LFImplementation<Physics>::compute_adv_residuals_no_flux_metric(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    interpolation::CellGeoMetric<MD1> const &CGM, interpolation::CellSolutionMetric<MD2> const &CSM,
    typename Physics::FluxV const &cell_residual, PGRDLFMethodData<Physics> &MD)
{
  // Reset the cell residuals
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n].fill(0.0);
    MD.m_LF_proj_vec[n].fill(0.0);
  }

  // Reset the cell update coefficient
  MD.m_k_LF.fill(0.0);
  MD.m_elem_wave_speed.fill(0.0);

  MD.m_u_avg.fill(0.0);

  // sol_nodal_values is a matrix of dimension
  // [(nb. nodes) x (nb. equations)]

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_u_avg += sol_nodal_values.row_transpose(n);
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    MD.m_u_avg[eq] *= 1. / MD.m_nb_nodes;
  }

  // Quadrature weights
  const math::DenseDVec<Real> &w = CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;
  // Solution in quadrature points

  const math::DenseConstMatView<Real> uq = CSM.field_values();

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = CSM.field_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // ==============================================================
  // PART I: COMPUTE THE (NONLIMITED) LF RESIDUALS (CENTRAL METHOD)
  // ==============================================================

  MD.m_cell_volume = 0.0;

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);
    // Compute the gradient of each shape function in physical space
    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    const Real wj_q = w[q] * jdet[q];
    MD.m_cell_volume += wj_q;

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      /*
      Physics::flux_jacobian_eigen_structure(MD.m_props,
      MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv, MD.Dvp);
      */

      Physics::flux_jacobian_eigen_values(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.m_eigenvalues);

      MD.m_max_eigenvalue = 1.e-9;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.m_eigenvalues[e]));
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));

        // THIS DOESN'T WORK (STABILITY PROBLEMS)!
        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      MD.m_k_LF[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

    } // Loop over nodes

  } // Loop over quadrature points

  // ----------------------------------------
  // Compute LF alpha (stability coefficient)
  // ----------------------------------------

  // At this point, MD.m_elem_wave_speed[n] holds the maximum spectral radius
  // of ( A \cdot \nabla \varphi_n ) sampled at quadrature points

  MD.m_LF_alpha = MD.m_k_LF[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_LF_alpha = std::max(std::abs(MD.m_k_LF[n]), MD.m_LF_alpha);
  }

  // -----------------------------------------
  // Finish computation of update coefficients
  // -----------------------------------------

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= MD.m_cell_volume;
  }

  // ---------------------------
  // Compute nodal residuals
  // ---------------------------
  MD.m_elem_node_res[0] = 1.0 / MD.m_nb_nodes * cell_residual;

  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] = MD.m_elem_node_res[0];
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n] += MD.m_LF_alpha * (sol_nodal_values.row_transpose(n) - MD.m_u_avg);
  }

  // Limiting for HO scheme
  residual_limiter_type::limit_cell_residuals(cell_residual, MD.m_sum_x_p, MD);
}

// ----------------------------------------------------------------------------

template <typename Physics>
struct CellSchemeSelector<Physics, LF>
{
  typedef LFImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
