#ifndef PDEKIT_Solver_RDM_PGB_hpp
#define PDEKIT_Solver_RDM_PGB_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "solver/rdm/RDMethodConstData.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/cellsplitters/CellSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Tag class to mark B scheme with variable beta coefficients

// ----------------------------------------------------------------------------

class PGB
{
};

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics>
class PGBImplementation
{
  public:
  enum
  {
    is_var_beta_type_rdm = 1
  };
  /// TYPEDEFS

  typedef Physics phys_model;
  typedef PGRDBMethodData<Physics> method_data;

  /// Constructor
  PGBImplementation();

  /// Destructor
  ~PGBImplementation();

  /// Compute the residuals on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      PGRDBMethodData<Physics> &MD);

  /// Compute the residuals and Jacobians on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_res_and_jacobian(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      PGRDBMethodData<Physics> &MD);

  /// Compute the residuals on one element, include source term
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_reaction_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      PGRDBMethodData<Physics> &MD);

  private:
  enum
  {
    NEQ = phys_model::NEQ
  };
};

// ----------------------------------------------------------------------------

template <typename Physics>
PGBImplementation<Physics>::PGBImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
PGBImplementation<Physics>::~PGBImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGBImplementation<Physics>::compute_adv_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    PGRDBMethodData<Physics> &MD)
{
  // Reset the element residuals
  MD.elem_node_res().fill(0.0);
  MD.n_node_res().fill(0.0);
  MD.lda_node_res().fill(0.0);

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

  Real max_nodal_theta = 0.0;

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    max_nodal_theta = std::max(max_nodal_theta, MD.m_blending_coeff[n]);
  }

  const bool skip_external_dissipation = MD.m_use_external_theta && (max_nodal_theta < 1.e-6);

#if 0

  /// ASSEMBLE THE NODAL RESIDUALS
  /// THIS IS A VERSION WHERE BLENDING COEFFICIENT IS INTEGRATED AT EACH QUADRATURE POINT
  /// IN THE END, WE TAKE THE MEAN VALUE OF THE BLENDING COEFF. IN CELL

  // Reset the blending coefficient
  MD.theta.fill(0.0);
  Real elem_volume = 0.0;

  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {

    MD.sum_Kp.fill(0.0);
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

    const math::ConstMatrixBlock<Real> inv_J = CGM.inv_jacobi(q);

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
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

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its absolute
        // value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.Dvp(e, e)));

        // MD.Dvm(e, e) = std::min(0.0, MD.Dvp(e, e));
        MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;
      // MD.m_Km[n] = MD.Rv * MD.Dvm * MD.Lv;
      MD.sum_Kp += MD.m_Kp[n];
    }

    MD.sum_Kp.inv(MD.inv_sum_Kp);


    /// a) Dissipation term
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      // MD.m_n_node_res[n].fill(0.0);

      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        if (m != n)
        {
          MD.m_n_node_res[n] +=
              wj_q * (MD.m_Kp[n] * MD.inv_sum_Kp * MD.m_Kp[m] *
                      (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
        }
      }
    }

    /// b) LDA scheme nodal residuals
    MD.m_res_at_point.fill(0.0);

    /*
    BOOST_FOREACH(const math::DynamicMatrix<Real> & dF, S.m_dF)
    {
      S.m_res_at_point += dF.row_transpose(q);
    }
    */

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // LDA part
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      // Here we DON'T ACCUMULATE, just store temporarily the LDA part of the
      // residual
      MD.m_lda_node_res[n] += wj_q * (MD.m_Kp[n] * MD.inv_sum_Kp * MD.m_res_at_point);
    }

    /// c) Compute the blending coefficient

    // MD.theta.fill(0.0);
    // Integrate the blending coefficient

    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      MD.abs_phi = 0.0;
      MD.abs_phiN = 0.0;

      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        MD.abs_phi += MD.m_lda_node_res[n][eq];
        MD.abs_phiN += std::abs(MD.m_lda_node_res[n][eq] + MD.m_n_node_res[n][eq]);
      }
      MD.theta(eq, eq) += wj_q * std::abs(MD.abs_phi) / (MD.abs_phiN + 1.e-6);
    }

    elem_volume += wj_q;

    //    /// d) Assemble LDA + theta * dissipation
    //    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    //    {
    //      for (Uint eq = 0; eq < NEQ; ++eq)
    //      {
    //        MD.m_elem_node_res[n][eq] +=
    //            wj_q * (MD.m_lda_node_res[n][eq] + MD.theta(eq, eq) * MD.m_n_node_res[n][eq]);
    //      }
    //    }
  } // Loop over quadrature points

  // Compute the mean value of blending coefficient
  for (Uint eq = 0; eq < NEQ; ++eq)
  {
    MD.theta(eq, eq) /= elem_volume;
  }

  /// d) Assemble LDA + theta * dissipation
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      MD.m_elem_node_res[n][eq] =
          (MD.m_lda_node_res[n][eq] + MD.theta(eq, eq) * MD.m_n_node_res[n][eq]);
    }
  }

  /// Update the wave speeds - multiply by element volume
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }

#endif

  // ------------------------

#if 1
  Real elem_volume = 0.0;

  /// ASSEMBLE THE NODAL RESIDUALS
  /// THIS IS A VERSION WHERE BLENDING COEFFICIENT IS COMPUTED CELL-WISE
  using JM = typename phys_model::JM;
  JM &Rv   = MD.Rv();
  JM &Lv   = MD.Lv();
  JM &Dvp  = MD.Dvp();
  // JM &Dvm = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();

  MD.m_flux_integral_in_elem.fill(0.0);

  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    sum_Kp.fill(0.0);
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

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, Rv, Lv, Dvp);

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        // MD.Dvm(e, e) = std::min(0.0, MD.Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      Kp[n] = Rv * Dvp * Lv;
      // MD.m_Km[n] = MD.Rv * MD.Dvm * MD.Lv;
      sum_Kp += Kp[n];
    } // Loop over the nodes of the element

    sum_Kp.inv(inv_sum_Kp);

    /// a) Dissipation term residuals
    if (!skip_external_dissipation)
    {
      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        for (Uint m = 0; m < MD.nb_nodes(); ++m)
        {
          if (m != n)
          {
            math::DenseVecView<Real> n_res_in_node = MD.n_node_res_block(n * NEQ, NEQ);
            n_res_in_node +=
                wj_q * (Kp[n] * inv_sum_Kp * Kp[m] *
                        (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
          }
        }
      }
    }

    /// b) LDA scheme nodal residuals
    MD.m_res_at_point.fill(0.0);

    /*
    BOOST_FOREACH(const math::DynamicMatrix<Real> & dF, S.m_dF)
    {
      S.m_res_at_point += dF.row_transpose(q);
    }
    */

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    MD.m_flux_integral_in_elem += wj_q * MD.m_res_at_point;

    // LDA part
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);

      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      // Here we ACCUMULATE THE LDA RESIDUAL
      lda_res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
    }

  } // Loop over quadrature points

  /// c) AFTER the numerical integration, compute blending coefficient

  if (!MD.m_use_external_theta)
  {
    MD.theta_cellwise.fill(0.0);
    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      MD.abs_phi  = 0.0;
      MD.abs_phiN = 0.0;

      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        // MD.abs_phi += MD.m_n_node_res[n][eq];
        // MD.abs_phiN += std::abs(MD.m_n_node_res[n][eq]);

        // MD.abs_phi += MD.m_lda_node_res[n][eq];
        math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(n * NEQ, NEQ);
        math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);

        MD.abs_phiN += std::abs(lda_res_in_node[eq] + n_res_in_node[eq]);
      }
      MD.theta_cellwise[eq] = std::abs(MD.m_flux_integral_in_elem[eq]) / (MD.abs_phiN + 1.e-6);
    }

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      MD.m_blending_coeff[n] = MD.theta_cellwise[0];
    }
  }

  /// Blend the dissipation term with the LDA
  /// nodal residuals
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(n * NEQ, NEQ);
    math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);
    math::DenseVecView<Real> b_res_in_node   = MD.elem_node_res_block(n * NEQ, NEQ);

    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      b_res_in_node[eq] = lda_res_in_node[eq] + MD.m_blending_coeff[n] * n_res_in_node[eq];
    }
  }

  /// Update the wave speeds - multiply by element volume
  /*
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.nb_nodes(); ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }

#endif
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGBImplementation<Physics>::compute_adv_res_and_jacobian(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    PGRDBMethodData<Physics> &MD)
{
  // Reset element jacobian
  MD.elem_jacobian().fill(0.0);
  MD.elem_jacobian_stab().fill(0.0);

  // Reset the element residuals
  MD.elem_node_res().fill(0.0);
  MD.n_node_res().fill(0.0);
  MD.lda_node_res().fill(0.0);

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

  Real max_nodal_theta = 0.0;

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    max_nodal_theta = std::max(max_nodal_theta, MD.m_blending_coeff[n]);
  }

  const bool skip_external_dissipation = MD.m_use_external_theta && (max_nodal_theta < 1.e-6);

  Real elem_volume = 0.0;
  typename Physics::JM beta_distr_mat;

  /// ASSEMBLE THE NODAL RESIDUALS
  /// THIS IS A VERSION WHERE BLENDING COEFFICIENT IS COMPUTED CELL-WISE
  using JM = typename phys_model::JM;
  JM &Rv   = MD.Rv();
  JM &Lv   = MD.Lv();
  JM &Dvp  = MD.Dvp();
  // JM &Dvm = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();

  MD.m_flux_integral_in_elem.fill(0.0);

  math::DenseMatView<Real> elem_jacobian      = MD.elem_jacobian();
  math::DenseMatView<Real> elem_jacobian_stab = MD.elem_jacobian_stab();

  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {
    sum_Kp.fill(0.0);
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

    for (Uint m = 0; m < MD.nb_nodes(); ++m)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, m);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, Rv, Lv, Dvp);

      // Here MD.m_Km[m] stores the matrix [K] corresponding to direction
      // given by grad(phi_m) This is needed to build the jacobians of the
      // LDA scheme
      Km[m] = Rv * Dvp * Lv;

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        // MD.Dvm(e, e) = std::min(0.0, MD.Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[m] = std::max(MD.m_elem_wave_speed[m], MD.m_max_eigenvalue);

      Kp[m] = Rv * Dvp * Lv;
      // MD.m_Km[n] = MD.Rv * MD.Dvm * MD.Lv;
      sum_Kp += Kp[m];
    } // Loop over m

    sum_Kp.inv(inv_sum_Kp);

    /// I) LDA scheme nodal residuals
    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    MD.m_flux_integral_in_elem += wj_q * MD.m_res_at_point;

    // LDA part
    for (Uint m = 0; m < MD.nb_nodes(); ++m)
    {
      // a) Accumulate into residuals (local RHS)
      math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(m * NEQ, NEQ);

      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      // Here we ACCUMULATE THE LDA RESIDUAL

      beta_distr_mat = Kp[m] * inv_sum_Kp;
      lda_res_in_node += wj_q * (beta_distr_mat * MD.m_res_at_point);

      // b) Accumulate into local Jacobi matrix (local LHS)
      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        // MD.m_jacobi_at_point = MD.m_Kp[m] * MD.inv_sum_Kp *
        // MD.m_Km[n];
        MD.m_jacobi_at_point = wj_q * (beta_distr_mat * Km[n]);

        for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
        {
          for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
          {
            elem_jacobian(m * NEQ + eq_m, n * NEQ + eq_n) += MD.m_jacobi_at_point(eq_m, eq_n);
          }
        }
      } // Loop over n
    }   // Loop over m

    /// II) Dissipation term residuals

    if (!skip_external_dissipation)
    {
      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        beta_distr_mat = Kp[m] * inv_sum_Kp;
        for (Uint n = 0; n < MD.nb_nodes(); ++n)
        {
          // Note that for m != n, the expression below is
          // Kplus(m) * inv(sum(Kplus)) * Kplus(n),
          // else it is
          // Kplus(m) * inv(sum(Kplus)) * Kplus(m)
          MD.m_jacobi_at_point = beta_distr_mat * Kp[n];

          if (m != n)
          {
            math::DenseVecView<Real> n_res_in_node = MD.n_node_res_block(m * NEQ, NEQ);
            n_res_in_node += wj_q * (MD.m_jacobi_at_point * (sol_nodal_values.row_transpose(m) -
                                                             sol_nodal_values.row_transpose(n)));

            for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
            {
              for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
              {
                elem_jacobian_stab(m * NEQ + eq_m, n * NEQ + eq_n) -=
                    wj_q * MD.m_jacobi_at_point(eq_m, eq_n);
              }
            }
          } // if ( m != n )
          else
          {
            // Contribution to rhs is zero for m == n
            // No contribution to rhs, only to Jacobian. It is:
            // Kplus(m) - { Kplus(m) * inv(sum(Kplus)) * Kplus(m) }
            for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
            {
              for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
              {
                elem_jacobian_stab(m * NEQ + eq_m, n * NEQ + eq_n) +=
                    wj_q * (Kp[m](eq_m, eq_n) - MD.m_jacobi_at_point(eq_m, eq_n));
              }
            }
          }
        } // Loop over n
      }   // Loop over m

    } // If the dissipation term has to be used

  } // Loop over quadrature points

  // ----------------------------------------------

  /// III) AFTER the numerical integration, compute blending coefficient

  // If the coefficient was computed externally, just fill theta
  if (!MD.m_use_external_theta)
  {
    MD.theta_cellwise.fill(0.0);
    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      MD.abs_phi  = 0.0;
      MD.abs_phiN = 0.0;

      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        // MD.abs_phi += MD.m_n_node_res[n][eq];
        // MD.abs_phiN += std::abs(MD.m_n_node_res[n][eq]);

        // MD.abs_phi += MD.m_lda_node_res[n][eq];
        math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(n * NEQ, NEQ);
        math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);

        MD.abs_phiN += std::abs(lda_res_in_node[eq] + n_res_in_node[eq]);
      }
      MD.theta_cellwise[eq] = std::abs(MD.m_flux_integral_in_elem[eq]) / (MD.abs_phiN + 1.e-6);
    }

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      MD.m_blending_coeff[n] = MD.theta_cellwise[0];
    }
  }

  // Debugging: switch off stabilization term completely
  // MD.theta.fill(0.0);

  /// Blend the dissipation term with the LDA
  /// nodal residuals
  for (Uint m = 0; m < MD.nb_nodes(); ++m)
  {
    math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(m * NEQ, NEQ);
    math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(m * NEQ, NEQ);
    math::DenseVecView<Real> b_res_in_node   = MD.elem_node_res_block(m * NEQ, NEQ);

    for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
    {
      b_res_in_node[eq_m] = lda_res_in_node[eq_m] + MD.m_blending_coeff[m] * n_res_in_node[eq_m];
    }
  }

  if (!skip_external_dissipation)
  {
    // Add Jacobian of stabilization to the LDA Jacobian
    for (Uint m = 0; m < MD.nb_nodes(); ++m)
    {
      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
        {
          for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
          {
            elem_jacobian(m * NEQ + eq_m, n * NEQ + eq_n) +=
                MD.m_blending_coeff[m] * elem_jacobian_stab(m * NEQ + eq_m, n * NEQ + eq_n);
          }
        }
      } // Loop over n
    }   // Loop over m

  } //  if !skip_external_dissipation

  /// Update the wave speeds - multiply by element volume
  /*
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.nb_nodes(); ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGBImplementation<Physics>::compute_adv_reaction_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    PGRDBMethodData<Physics> &MD)
{
  // Reset the element residuals
  MD.elem_node_res().fill(0.0);
  MD.n_node_res().fill(0.0);
  MD.lda_node_res().fill(0.0);

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

  Real max_nodal_theta = 0.0;

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    max_nodal_theta = std::max(max_nodal_theta, MD.m_blending_coeff[n]);
  }

  const bool skip_external_dissipation = MD.m_use_external_theta && (max_nodal_theta < 1.e-6);

  Real elem_volume = 0.0;

  /// ASSEMBLE THE NODAL RESIDUALS
  /// THIS IS A VERSION WHERE BLENDING COEFFICIENT IS COMPUTED CELL-WISE
  using JM = typename phys_model::JM;
  JM &Rv   = MD.Rv();
  JM &Lv   = MD.Lv();
  JM &Dvp  = MD.Dvp();
  // JM &Dvm = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  // common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();

  MD.m_flux_integral_in_elem.fill(0.0);

  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    sum_Kp.fill(0.0);
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

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, Rv, Lv, Dvp);

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        // MD.Dvm(e, e) = std::min(0.0, MD.Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      Kp[n] = Rv * Dvp * Lv;
      // MD.m_Km[n] = MD.Rv * MD.Dvm * MD.Lv;
      sum_Kp += Kp[n];
    } // Loop over the nodes of the element

    sum_Kp.inv(inv_sum_Kp);

    /// a) Dissipation term residuals
    if (!skip_external_dissipation)
    {
      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        for (Uint m = 0; m < MD.nb_nodes(); ++m)
        {
          if (m != n)
          {
            math::DenseVecView<Real> n_res_in_node = MD.n_node_res_block(n * NEQ, NEQ);

            n_res_in_node +=
                wj_q * (Kp[n] * inv_sum_Kp * Kp[m] *
                        (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
          }
        }
      }
    }

    /// b) LDA scheme nodal residuals
    MD.m_res_at_point.fill(0.0);

    /*
    BOOST_FOREACH(const math::DynamicMatrix<Real> & dF, S.m_dF)
    {
      S.m_res_at_point += dF.row_transpose(q);
    }
    */

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    // Include source term
    MD.m_res_at_point -= sq.row_transpose(q);

    MD.m_flux_integral_in_elem += wj_q * MD.m_res_at_point;

    // LDA part
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);

      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      // Here we ACCUMULATE THE LDA RESIDUAL
      lda_res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
    }

  } // Loop over quadrature points

  /// c) AFTER the numerical integration, compute blending coefficient
  if (!MD.m_use_external_theta)
  {
    MD.theta_cellwise.fill(0.0);
    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      MD.abs_phi  = 0.0;
      MD.abs_phiN = 0.0;

      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        // MD.abs_phi += MD.m_n_node_res[n][eq];
        // MD.abs_phiN += std::abs(MD.m_n_node_res[n][eq]);

        // MD.abs_phi += MD.m_lda_node_res[n][eq];
        math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(n * NEQ, NEQ);
        math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);

        MD.abs_phiN += std::abs(lda_res_in_node[eq] + n_res_in_node[eq]);
      }
      MD.theta_cellwise[eq] = std::abs(MD.m_flux_integral_in_elem[eq]) / (MD.abs_phiN + 1.e-6);
    }

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      MD.m_blending_coeff[n] = MD.theta_cellwise[0];
    }
  }

  /// Blend the dissipation term with the LDA
  /// nodal residuals
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    math::DenseVecView<Real> n_res_in_node   = MD.n_node_res_block(n * NEQ, NEQ);
    math::DenseVecView<Real> lda_res_in_node = MD.lda_node_res_block(n * NEQ, NEQ);
    math::DenseVecView<Real> b_res_in_node   = MD.elem_node_res_block(n * NEQ, NEQ);

    for (Uint eq = 0; eq < NEQ; ++eq)
    {
      b_res_in_node[eq] = lda_res_in_node[eq] + MD.m_blending_coeff[n] * n_res_in_node[eq];
    }
  }

  /// Update the wave speeds - multiply by element volume
  /*
  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] *= elem_volume;
  }
  */

  Real max_wave_speed = MD.m_elem_wave_speed[0];
  for (Uint n = 1; n < MD.nb_nodes(); ++n)
  {
    max_wave_speed = std::max(max_wave_speed, MD.m_elem_wave_speed[n]);
  }

  for (Uint n = 0; n < MD.nb_nodes(); ++n)
  {
    MD.m_elem_wave_speed[n] = elem_volume * max_wave_speed;
  }
}

// ----------------------------------------------------------------------------

template <typename Physics>
struct CellSchemeSelector<Physics, PGB>
{
  typedef PGBImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
