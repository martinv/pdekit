#ifndef PDEKIT_Solver_RDM_PGLF1_hpp
#define PDEKIT_Solver_RDM_PGLF1_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/cellsplitters/LFResidualLimiter.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename Physics>
class PGLF1
{
  public:
  /// TYPEDEFS

  typedef Physics phys_model;
  typedef PGRDLFMethodData<Physics> method_data;

  /// Constructor
  PGLF1();

  /// Destructor
  ~PGLF1();

  /// Compute the residuals on one element
  template <typename MD1, typename MD2, typename MD3>
  static void compute_adv_residuals(math::DenseConstMatView<Real> const &sol_nodal_values,
                                    interpolation::CellGeoMetric<MD1> const &CGM,
                                    interpolation::CellSolutionMetric<MD2> const &CSM,
                                    interpolation::CellFluxMetric<MD3, Physics> const &CFM,
                                    PGRDLFMethodData<Physics> &MD);

  /// Compute the residuals on one element, include source term
  template <typename MD1, typename MD2, typename MD3>
  static void compute_adv_reaction_residuals(math::DenseConstMatView<Real> const &sol_nodal_values,
                                             interpolation::CellGeoMetric<MD1> const &CGM,
                                             interpolation::CellSolutionMetric<MD2> const &CSM,
                                             interpolation::CellFluxMetric<MD3, Physics> const &CFM,
                                             interpolation::CellSolutionMetric<MD2> const &CSrcM,
                                             PGRDLFMethodData<Physics> &MD);

  private:
  /// TYPES
  typedef typename internal::LFResidualLimiter<Physics, Physics::NEQ> residual_limiter_type;
};

// ----------------------------------------------------------------------------

template <typename Physics>
PGLF1<Physics>::PGLF1()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
PGLF1<Physics>::~PGLF1()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3>
void PGLF1<Physics>::compute_adv_residuals(math::DenseConstMatView<Real> const &sol_nodal_values,
                                           interpolation::CellGeoMetric<MD1> const &CGM,
                                           interpolation::CellSolutionMetric<MD2> const &CSM,
                                           interpolation::CellFluxMetric<MD3, Physics> const &CFM,
                                           PGRDLFMethodData<Physics> &MD)
{

  // Reset the element residuals
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n].fill(0.0);
  }

  // Reset the element update coefficient
  MD.m_elem_wave_speed.fill(0.0);

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
  // std::cout << "Solution at quadrature points = " << std::endl << uq <<
  // std::endl;

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = CSM.field_derivatives(d);
    MD.m_grad_F[d] = CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // std::cout << std::endl << std::endl;

  const Real one_over_nb_nodes = 1.0 / MD.m_nb_nodes;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    const Real wj_q = w[q] * jdet[q];

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);

    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    MD.m_LF_alpha = 1.e-9;

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv,
                                             MD.Dvp);

      MD.m_max_eigenvalue = 1.e-9;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.Dvp(e, e)));
      }

      MD.m_elem_wave_speed[n] += wj_q * (one_over_nb_nodes + MD.m_max_eigenvalue);

      MD.m_LF_alpha = std::max(MD.m_LF_alpha, MD.m_max_eigenvalue);

    } // Loop over nodes of the element

    MD.m_res_at_point.fill(0.0);

    /*
    BOOST_FOREACH(const math::DynamicMatrix<Real> & dF, MD.m_dF)
    {
      MD.m_res_at_point += dF.row_transpose(q);
    }
    */

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Central part
    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      // vector = scalar *  scalar * vector
      // [NEQ]  = scalar *  scalar * [NEQ]
      MD.m_elem_node_res[n] += wj_q * MD.m_res_at_point;
    }

    // Stabilization

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      for (Uint m = 0; m < MD.m_nb_nodes; ++m)
      {
        if (m != n)
        {
          MD.m_elem_node_res[n] += (wj_q * MD.m_LF_alpha) * (sol_nodal_values.row_transpose(n) -
                                                             sol_nodal_values.row_transpose(m));
        }
      }
    }

  } // Loop over quadrature points

  // Limiting - works only in the scalar case for the moment
  MD.m_total_cell_res.fill(0.0);

  // Compute the total residual in the cell
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    // Divide the nodal residual by number of local nodes in element
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      MD.m_elem_node_res[n][eq] *= one_over_nb_nodes;
    }
    MD.m_total_cell_res += MD.m_elem_node_res[n];
  }

  // Limiting for HO scheme
  residual_limiter_type::limit_cell_residuals(MD.m_total_cell_res, MD.m_sum_x_p, MD);
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3>
void PGLF1<Physics>::compute_adv_reaction_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    interpolation::CellGeoMetric<MD1> const &CGM, interpolation::CellSolutionMetric<MD2> const &CSM,
    interpolation::CellFluxMetric<MD3, Physics> const &CFM,
    interpolation::CellSolutionMetric<MD2> const &CSrcM, PGRDLFMethodData<Physics> &MD)
{
  // Reset the element residuals
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_node_res[n].fill(0.0);
  }

  // Reset the element update coefficient
  MD.m_elem_wave_speed.fill(0.0);

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
  // std::cout << "Solution at quadrature points = " << std::endl << uq <<
  // std::endl;

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

  // std::cout << std::endl << std::endl;

  const Real one_over_nb_nodes = 1.0 / MD.m_nb_nodes;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    const Real wj_q = w[q] * jdet[q];

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);

    const math::DenseConstMatView<Real> inv_J = CGM.inv_jacobi(q);

    MD.m_LF_alpha = 1.e-9;

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv,
                                             MD.Dvp);

      MD.m_max_eigenvalue = 1.e-9;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.Dvp(e, e)));
      }

      MD.m_elem_wave_speed[n] += wj_q * (one_over_nb_nodes + MD.m_max_eigenvalue);

      MD.m_LF_alpha = std::max(MD.m_LF_alpha, MD.m_max_eigenvalue);

    } // Loop over the nodes of the element

    MD.m_res_at_point.fill(0.0);

    /*
    BOOST_FOREACH(const math::DynamicMatrix<Real> & dF, MD.m_dF)
    {
      MD.m_res_at_point += dF.row_transpose(q);
    }
    */

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Include source term
    MD.m_res_at_point -= sq.row_transpose(q);

    // Central part
    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      // vector = scalar *  scalar * vector
      // [NEQ]  = scalar *  scalar * [NEQ]
      MD.m_elem_node_res[n] += wj_q * (wj_q * one_over_nb_nodes) * MD.m_res_at_point;
    }

    // Stabilization

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      for (Uint m = 0; m < MD.m_nb_nodes; ++m)
      {
        if (m != n)
        {
          MD.m_elem_node_res[n] +=
              (wj_q * one_over_nb_nodes * MD.m_LF_alpha) *
              (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m));
        }
      }
    }

  } // Loop over quadrature points

  // Limiting - works only in the scalar case for the moment
  MD.m_total_cell_res.fill(0.0);

  // Compute the total residual in the cell
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_total_cell_res += MD.m_elem_node_res[n];
  }

  MD.m_sum_x_p.fill(0.0);
  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    if (MD.m_total_cell_res[eq] != 0.0)
    {
      for (Uint n = 0; n < MD.m_nb_nodes; ++n)
      {
        // Compute x+ and store it (temporarily) in
        // MD.m_elemnode_res[n][eq]
        MD.m_elem_node_res[n][eq] =
            std::max(MD.m_elem_node_res[n][eq] / MD.m_total_cell_res[eq], 0.0);
        MD.m_sum_x_p[eq] += MD.m_elem_node_res[n][eq];
      }

      for (Uint n = 0; n < MD.m_nb_nodes; ++n)
      {
        MD.m_elem_node_res[n][eq] =
            MD.m_elem_node_res[n][eq] / MD.m_sum_x_p[eq] * MD.m_total_cell_res[eq];
      }

    } // If total cell residual for this equation != 0
    else
    {
      for (Uint n = 0; n < MD.m_nb_nodes; ++n)
      {
        MD.m_elem_node_res[n][eq] = 0.0;
      }
    }
  }
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
