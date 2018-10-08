#ifndef PDEKIT_Solver_RDM_PGN_hpp
#define PDEKIT_Solver_RDM_PGN_hpp

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

// Tag class to mark PGN with variable beta coefficients

// ----------------------------------------------------------------------------

class PGN
{
};

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics>
class PGNImplementation
{
  public:
  enum
  {
    is_var_beta_type_rdm = 1
  };
  /// TYPEDEFS

  typedef Physics phys_model;
  typedef RDMethodScratchData<Physics> method_data;

  /// Constructor
  PGNImplementation();

  /// Destructor
  ~PGNImplementation();

  /// Compute the residuals on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      RDMethodScratchData<Physics> &MD);

  /// Compute the residuals and Jacobians on one element
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_res_and_jacobian(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      RDMethodScratchData<Physics> &MD);

  /// Compute the residuals on one element, include source term
  template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
  static void compute_adv_reaction_residuals(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
      RDMethodScratchData<Physics> &MD);

  private:
  enum
  {
    NEQ = phys_model::NEQ
  };
};

// ----------------------------------------------------------------------------

template <typename Physics>
PGNImplementation<Physics>::PGNImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
PGNImplementation<Physics>::~PGNImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGNImplementation<Physics>::compute_adv_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    RDMethodScratchData<Physics> &MD)
{
  // Reset the element residuals
  MD.elem_node_res().fill(0.0);

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

  /// ASSEMBLE THE NODAL RESIDUALS
  ///   typename phys_model::JM &Rv = MD.Rv();
  using JM       = typename phys_model::JM;
  JM &Rv         = MD.Rv();
  JM &Lv         = MD.Lv();
  JM &Dvp        = MD.Dvp();
  JM &Dvm        = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();

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

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        Dvm(e, e) = std::min(0.0, Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      Kp[n] = Rv * Dvp * Lv;
      Km[n] = Rv * Dvm * Lv;
      sum_Kp += Kp[n];
    }

    sum_Kp.inv(inv_sum_Kp);

    //#define Plain_N

#ifdef Plain_N
    /// a) Just N scheme
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        if (m != n)
        {
          MD.m_elem_node_res[n] -=
              wj_q * (MD.m_Kp[n] * MD.inv_sum_Kp * MD.m_Km[m] *
                      (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
        }
      }
    }
#endif

#define LDA_and_dissipation

#ifdef LDA_and_dissipation
    /// b) LDA + dissipation
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

    // LDA part
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
    }

    // Dissipation

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        if (m != n)
        {
          math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
          res_in_node +=
              wj_q * (Kp[n] * inv_sum_Kp * Kp[m] *
                      (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
        }
      }
    }
#endif
  } // Loop over quadrature points

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
void PGNImplementation<Physics>::compute_adv_res_and_jacobian(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    RDMethodScratchData<Physics> &MD)
{
  // Reset element jacobian
  MD.elem_jacobian().fill(0.0);

  // Reset the element residuals
  MD.elem_node_res().fill(0.0);

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
  typename Physics::JM K_dot_sumK_inv;

  math::DenseMatView<Real> elem_jacobian = MD.elem_jacobian();

  /// ASSEMBLE THE NODAL RESIDUALS
  using JM       = typename phys_model::JM;
  JM &Rv         = MD.Rv();
  JM &Lv         = MD.Lv();
  JM &Dvp        = MD.Dvp();
  JM &Dvm        = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();

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

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        Dvm(e, e) = std::min(0.0, Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      Kp[n] = Rv * Dvp * Lv;
      Km[n] = Rv * Dvm * Lv;
      // Note that sum_Kp holds the sum of K^{-} here!
      sum_Kp += Km[n];
    }

    sum_Kp.inv(inv_sum_Kp);

    /// Just N scheme
    for (Uint m = 0; m < MD.nb_nodes(); ++m)
    {
      // Temporary K_m^{+} * sum(K^{-})
      K_dot_sumK_inv = Kp[m] * inv_sum_Kp;
      for (Uint n = 0; n < MD.nb_nodes(); ++n)
      {
        if (m != n)
        {
          math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(m * NEQ, NEQ);
          res_in_node +=
              wj_q * (Kp[m] * inv_sum_Kp * Km[n] *
                      (sol_nodal_values.row_transpose(m) - sol_nodal_values.row_transpose(n)));

          MD.m_jacobi_at_point = -1. * (K_dot_sumK_inv * Km[n]);
        }
        else
        {
          MD.m_jacobi_at_point = Kp[m] - K_dot_sumK_inv * Km[n];
        }

        for (Uint eq_m = 0; eq_m < NEQ; ++eq_m)
        {
          for (Uint eq_n = 0; eq_n < NEQ; ++eq_n)
          {
            elem_jacobian(m * NEQ + eq_m, n * NEQ + eq_n) +=
                wj_q * MD.m_jacobi_at_point(eq_m, eq_n);
          }
        }
      } // Loop over n
    }   // Loop over m
  }     // Loop over quadrature points

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
void PGNImplementation<Physics>::compute_adv_reaction_residuals(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    const RDMethodConstData<CellGeoMetType, CellSolMetType, CellFluxMetType> &ConstMD,
    RDMethodScratchData<Physics> &MD)
{
  // Reset the element residuals
  MD.elem_node_res().fill(0.0);

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

  /// ASSEMBLE THE NODAL RESIDUALS
  using JM       = typename phys_model::JM;
  JM &Rv         = MD.Rv();
  JM &Lv         = MD.Lv();
  JM &Dvp        = MD.Dvp();
  JM &Dvm        = MD.Dvm();
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();

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

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // Here, MD.Dvp(e,e) is not positive yet, so we have to take its
        // absolute value to get spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(Dvp(e, e)));

        Dvm(e, e) = std::min(0.0, Dvp(e, e));
        Dvp(e, e) = std::max(0.0, Dvp(e, e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      Kp[n] = Rv * Dvp * Lv;
      Km[n] = Rv * Dvm * Lv;
      sum_Kp += Kp[n];
    }

    sum_Kp.inv(inv_sum_Kp);

    //#define Plain_N

#ifdef Plain_N
    /// FIXME: ADD SOURCE TERM HERE!
    /// a) Just N scheme
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        if (m != n)
        {
          MD.m_elem_node_res[n] -=
              wj_q * (MD.m_Kp[n] * MD.inv_sum_Kp * MD.m_Km[m] *
                      (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
        }
      }
    }
#endif

#define LDA_and_dissipation

#ifdef LDA_and_dissipation
    /// b) LDA + dissipation
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

    // LDA part
    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
      // vector = scalar *  (   matrix   *  matrix   * vector )
      // [NEQ]  = scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
    }

    // Dissipation

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      for (Uint m = 0; m < MD.nb_nodes(); ++m)
      {
        if (m != n)
        {
          math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
          res_in_node +=
              wj_q * (Kp[n] * inv_sum_Kp * Kp[m] *
                      (sol_nodal_values.row_transpose(n) - sol_nodal_values.row_transpose(m)));
        }
      }
    }
#endif
  } // Loop over quadrature points

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
struct CellSchemeSelector<Physics, PGN>
{
  typedef PGNImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
