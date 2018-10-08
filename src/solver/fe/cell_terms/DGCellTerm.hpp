#ifndef PDEKIT_Solver_DG_Cell_Term_hpp
#define PDEKIT_Solver_DG_Cell_Term_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "solver/fe/DGMethodConstData.hpp"
#include "solver/fe/DGMethodScratchData.hpp"

namespace pdekit
{

namespace solver
{

namespace fe
{

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics>
class DGCellTerm
{
  public:
  /// TYPEDEFS

  typedef Physics phys_model;
  typedef DGMethodScratchData<Physics> method_data;

  /// Constructor
  DGCellTerm();

  /// Destructor
  ~DGCellTerm();

  /// Compute the residuals on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      DGMethodScratchData<Physics> &MD);

  /// Compute the residuals and Jacobians on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_res_and_jacobian(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      DGMethodScratchData<Physics> &MD);

  /// Compute the residuals on one element, include source term
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_reaction_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      DGMethodScratchData<Physics> &MD);

  private:
  enum
  {
    NEQ = phys_model::NEQ
  };
};

// ----------------------------------------------------------------------------

template <typename Physics>
DGCellTerm<Physics>::DGCellTerm()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
DGCellTerm<Physics>::~DGCellTerm()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void DGCellTerm<Physics>::compute_adv_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    DGMethodScratchData<Physics> &MD)
{

  // Reset the element residuals
  MD.m_elem_node_res.fill(0.0);

  // Reset the element update coefficient
  MD.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &w = ConstMD.CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = ConstMD.CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = ConstMD.CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Solution in quadrature points
  const math::DenseConstMatView<Real> uq = ConstMD.CSM.field_values();
  // std::cout << "Solution at quadrature points = " << std::endl << uq <<
  // std::endl;

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = ConstMD.CSM.field_derivatives(d);
    MD.m_grad_F[d] = ConstMD.CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // std::cout << std::endl << std::endl;

  Real elem_volume = 0.0;

  typename Physics::FluxV eigenvalues;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    MD.sum_Kp.fill(0.0);
    const Real wj_q = w[q] * jdet[q];

    elem_volume += wj_q;

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);

    const math::DenseConstMatView<Real> inv_J = ConstMD.CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_values(MD.m_props, MD.m_grad_sf_at_pt_phys, eigenvalues);

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(eigenvalues[e]));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      /*
      MD.m_elem_wave_speed[n] += wj_q * (1. / MD.m_nb_nodes +
      MD.m_max_eigenvalue); MD.m_elem_wave_speed[n] = std::max(1. /
      MD.m_nb_nodes + MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);
      */

      MD.m_elem_wave_speed[n] =
          std::max(MD.m_elem_wave_speed[n], 1. / MD.m_nb_nodes + MD.m_max_eigenvalue);
    }

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

    // Diagonal matrix of shape function values
    // of shape function n at quadrature point q

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      MD.m_sf_at_pt.fill(0.0);
      for (Uint eq = 0; eq < NEQ; ++eq)
      {
        MD.m_sf_at_pt(eq, eq) = MD.m_V(q, n);
      }

      math::DenseVecView<Real> res_in_node = MD.m_elem_node_res.block(n * NEQ, NEQ);

      // vector   += scalar * (   matrix   *  matrix   * vector )
      // [NEQ]    += scalar * ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (MD.m_sf_at_pt * MD.m_res_at_point);
      //      S.m_elem_node_res[n] += wj_q * ((S.m_sf_at_pt + S.inv_sum_Kp
      //      *
      // S.m_Kp[n]) * S.m_res_at_point);
    }

  } // Loop over quadrature points

  /*
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void DGCellTerm<Physics>::compute_adv_res_and_jacobian(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    DGMethodScratchData<Physics> &MD)
{
  // Reset element jacobian
  MD.m_elem_jacobian.fill(0.0);

  // Reset the element residuals
  MD.m_elem_node_res.fill(0.0);

  // Reset the element update coefficient
  MD.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &w = ConstMD.CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = ConstMD.CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = ConstMD.CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Solution in quadrature points
  const math::DenseConstMatView<Real> uq = ConstMD.CSM.field_values();
  // std::cout << "Solution at quadrature points = " << std::endl << uq <<
  // std::endl;

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = ConstMD.CSM.field_derivatives(d);
    MD.m_grad_F[d] = ConstMD.CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // std::cout << std::endl << std::endl;

  Real elem_volume = 0.0;
  typename Physics::JM beta_distr_mat;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    MD.sum_Kp.fill(0.0);
    const Real wj_q = w[q] * jdet[q];

    elem_volume += wj_q;

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);

    const math::DenseConstMatView<Real> inv_J = ConstMD.CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv,
                                             MD.Dvp);

      // In SUPG, S.m_Kp stores the matrix K, not K+ !
      MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.Dvp(e, e)));
        MD.Dvp(e, e)        = std::max(0.0, MD.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      /*
      MD.m_elem_wave_speed[n] += wj_q * (1. / MD.m_nb_nodes +
      MD.m_max_eigenvalue); MD.m_elem_wave_speed[n] = std::max(1. /
      MD.m_nb_nodes + MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);
      */

      MD.m_elem_wave_speed[n] =
          std::max(MD.m_elem_wave_speed[n], 1. / MD.m_nb_nodes + MD.m_max_eigenvalue);

      MD.sum_Kp += (MD.Rv * MD.Dvp * MD.Lv);
    }

    MD.sum_Kp.inv(MD.inv_sum_Kp);

    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Diagonal matrix of shape function values
    // of shape function n at quadrature point q

    for (Uint m = 0; m < MD.m_nb_nodes; ++m)
    {
      // I) Accumulate into residuals (local RHS)

      beta_distr_mat.fill(0.0);
      for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
      {
        beta_distr_mat(eq_m, eq_m) = MD.m_V(q, m);
      }

      beta_distr_mat += MD.m_Kp[m] * MD.inv_sum_Kp;

      math::DenseVecView<Real> res_in_node = MD.m_elem_node_res.block(m * NEQ, NEQ);

      // vector   += scalar * (   matrix   *  matrix   * vector )
      // [NEQ]    += scalar * ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (beta_distr_mat * MD.m_res_at_point);
      //      S.m_elem_node_res[n] += wj_q * ((S.m_sf_at_pt + S.inv_sum_Kp
      //      *
      //                                       S.m_Kp[n]) *
      //                                       S.m_res_at_point);

      // II) Accumulate into local Jacobi matrix (local LHS)
      for (Uint n = 0; n < MD.m_nb_nodes; ++n)
      {
        // MD.m_jacobi_at_point = MD.m_Kp[m] * MD.inv_sum_Kp *
        // MD.m_Km[n];
        MD.m_jacobi_at_point = wj_q * (beta_distr_mat * MD.m_Kp[n]);

        for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
        {
          for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
          {
            MD.m_elem_jacobian(m * NEQ + eq_m, n * NEQ + eq_n) += MD.m_jacobi_at_point(eq_m, eq_n);
          }
        }
      } // Loop over n

    } // Loop over m

  } // Loop over quadrature points

  /*
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void DGCellTerm<Physics>::compute_adv_reaction_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const DGMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    DGMethodScratchData<Physics> &MD)
{
  // Reset the element residuals
  MD.m_elem_node_res.fill(0.0);

  // Reset the element update coefficient
  MD.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &w = ConstMD.CGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = ConstMD.CGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = ConstMD.CGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Solution in quadrature points
  const math::DenseConstMatView<Real> uq = ConstMD.CSM.field_values();
  // std::cout << "Solution at quadrature points = " << std::endl << uq <<
  // std::endl;

  // Source field in quadrature points
  const math::DenseConstMatView<Real> sq = ConstMD.CSrcM.field_values();

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MD.m_grad_u[d] = ConstMD.CSM.field_derivatives(d);
    MD.m_grad_F[d] = ConstMD.CFM.flux_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  // std::cout << std::endl << std::endl;

  Real elem_volume = 0.0;

  typename Physics::FluxV eigenvalues;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    MD.sum_Kp.fill(0.0);
    const Real wj_q = w[q] * jdet[q];

    elem_volume += wj_q;

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MD.m_grad_u_at_point.insert_col(dim, MD.m_grad_u[dim].row(q));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq.row_transpose(q), MD.m_grad_u_at_point,
                                MD.m_props);

    const math::DenseConstMatView<Real> inv_J = ConstMD.CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_values(MD.m_props, MD.m_grad_sf_at_pt_phys, eigenvalues);

      // In SUPG, S.m_Kp stores the matrix K, not K+ !
      MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(eigenvalues[e]));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (1. / MD.m_nb_nodes +
      // MD.m_max_eigenvalue);
      /*
      MD.m_elem_wave_speed[n] =
          std::max(1. / MD.m_nb_nodes + MD.m_elem_wave_speed[n],
      MD.m_max_eigenvalue);
      */

      MD.m_elem_wave_speed[n] =
          std::max(MD.m_elem_wave_speed[n], 1. / MD.m_nb_nodes + MD.m_max_eigenvalue);
    }

    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Option A) for source term: add it into the
    // residual - this whole residual will then be weighted
    // by the 'upwind' trial function
    MD.m_res_at_point -= sq.row_transpose(q);

    // Diagonal matrix of shape function values
    // of shape function n at quadrature point q

    for (Uint n = 0; n < MD.m_nb_nodes; ++n)
    {
      MD.m_sf_at_pt.fill(0.0);
      for (Uint eq = 0; eq < NEQ; ++eq)
      {
        MD.m_sf_at_pt(eq, eq) = MD.m_V(q, n);
      }

      math::DenseVecView<Real> res_in_node = MD.m_elem_node_res.block(n * NEQ, NEQ);

      // vector   += scalar * (   matrix   *  matrix   * vector )
      // [NEQ]    += scalar * ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (MD.m_sf_at_pt * MD.m_res_at_point);
      //      S.m_elem_node_res[n] += wj_q * ((S.m_sf_at_pt + S.inv_sum_Kp
      //      *
      // S.m_Kp[n]) * S.m_res_at_point);

      // Option B) for source term: use only Galerkin test function for
      // weighting of the source term MD.m_elem_node_res[n] -= wj_q *
      // MD.m_V(q,n) * sq.row_transpose(q);
    }

  } // Loop over quadrature points

  /*
  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.m_nb_nodes; ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.m_nb_nodes; ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }
}

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace fe

} // namespace solver

} // namespace pdekit

#endif
