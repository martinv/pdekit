#ifndef PDEKIT_Solver_RDM_PGLDA_hpp
#define PDEKIT_Solver_RDM_PGLDA_hpp

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

// Tag class to mark LDA with variable beta coefficients

// ----------------------------------------------------------------------------

class PGLDA
{
};

// ----------------------------------------------------------------------------

namespace internal
{

#define USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA 1

template <typename Physics>
class PGLDAImplementation
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
  PGLDAImplementation();

  /// Destructor
  ~PGLDAImplementation();

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

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2>
  static void compute_adv_residuals_no_flux_metric(
      math::DenseConstMatView<Real> const &sol_nodal_values,
      interpolation::CellGeoMetric<MD1> const &CGM,
      interpolation::CellSolutionMetric<MD2> const &CSM,
      typename Physics::FluxV const &cell_residual, RDMethodScratchData<Physics> &MD);

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
PGLDAImplementation<Physics>::PGLDAImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
PGLDAImplementation<Physics>::~PGLDAImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGLDAImplementation<Physics>::compute_adv_residuals(
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

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA
  math::DenseSVec<Real, NEQ> eigvalues;
#endif

  /// ASSEMBLE THE NODAL RESIDUALS
  using JM = typename phys_model::JM;
  /*
  JM &Rv = MD.Rv();
  JM &Lv = MD.Lv();
  JM &Dvp = MD.Dvp();
  JM &Dvm = MD.Dvm();
  */
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();

  for (Uint q = 0; q < ConstMD.CGM.nb_qd_pts(); ++q)
  {

    MD.sum_Kp().fill(0.0);
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

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA
      Physics::build_K_mat(MD.m_props, MD.m_grad_sf_at_pt_phys, eigvalues, Kp[n],
                           [](const Real eigvalue) { return (std::max(0.0, eigvalue)); });

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(eigvalues[e]));
      }

      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

#else
      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv,
                                             MD.Dvp);

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(MD.Dvp(e, e)));
        Dvp(e, e)           = std::max(0.0, Dvp(e, e));

        // Option b: largest positive eigenvalue
        // Dvp(e, e) = std::max(0.0, Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, Dvp(e,
        // e));
      }

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      MD.m_Kp[n] = MD.Rv * MD.Dvp * MD.Lv;
#endif

      sum_Kp += Kp[n];

      // ---------------------------------------------
      // Artificial viscosity term for shock capturing
      // ---------------------------------------------
      if (MD.m_art_visc[q] != 0.0)
      {
        math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
        for (Uint dim = 0; dim < Physics::DIM; ++dim)
        {
          res_in_node += (wj_q * MD.m_art_visc[q] * MD.m_grad_sf_at_pt_phys[dim]) *
                         MD.m_grad_u[dim].row_transpose(q);
        }
      }
    } // Loop over nodes of the element

    sum_Kp.inv(inv_sum_Kp);

    MD.m_res_at_point.fill(0.0);

    // If the physics has source term, compute it as a part of the residual
    phys_model::source_term(MD.m_props, MD.m_res_at_point);

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

    /*
    std::array<typename Physics::JM, Physics::DIM> flux_jacob;
    Physics::residual(MD.m_props, flux_jacob, MD.m_res_at_point);
    */

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
    }

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
void PGLDAImplementation<Physics>::compute_adv_res_and_jacobian(
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
    MD.m_grad_u[d] = ConstMD.CSM.field_derivatives(d);
    MD.m_grad_F[d] = ConstMD.CFM.flux_derivatives(d);
  }

  Real elem_volume = 0.0;
  typename Physics::JM beta_distr_mat;

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA
  math::DenseSVec<Real, NEQ> eigvalues;
#endif

  /// ASSEMBLE THE NODAL RESIDUALS
  using JM = typename phys_model::JM;
  /*
  JM &Rv = MD.Rv();
  JM &Lv = MD.Lv();
  JM &Dvp = MD.Dvp();
  JM &Dvm = MD.Dvm();
  */
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();
  common::ArrayView<typename phys_model::JM, _1D, Uint> Km = MD.Km();
  math::DenseMatView<Real> elem_jacobian                   = MD.elem_jacobian();

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

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA

      Physics::build_K_mat(MD.m_props, MD.m_grad_sf_at_pt_phys, eigvalues, Kp[n],
                           [](const Real eigvalue) { return (std::max(0.0, eigvalue)); });

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(eigvalues[e]));
      }

      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      // Here MD.m_Km[n] stores the matrix [K] corresponding to direction
      // given by grad(phi_n)
      Physics::build_K_mat(MD.m_props, MD.m_grad_sf_at_pt_phys, eigvalues, Km[n],
                           [](const Real eigvalue) { return (eigvalue); });

#else
      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, MD.Rv, MD.Lv,
                                             MD.Dvp);

      // Here MD.m_Km[n] stores the matrix [K] corresponding to direction
      // given by grad(phi_n)
      MD.m_Km[n] = MD.Rv * MD.Dvp * MD.Lv;

      // Now modify eigenvalues so that MD.Dvp really stores __positive__
      // eigenvalues
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

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      m_Kp[n] = Rv * Dvp * Lv;

#endif
      sum_Kp += Kp[n];

      // ---------------------------------------------
      // Artificial viscosity term for shock capturing
      // ---------------------------------------------
      if (MD.m_art_visc[q] != 0.0)
      {
        math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
        for (Uint dim = 0; dim < Physics::DIM; ++dim)
        {
          res_in_node += (wj_q * MD.m_art_visc[q] * MD.m_grad_sf_at_pt_phys[dim]) *
                         MD.m_grad_u[dim].row_transpose(q);
        }
      }
    } // Loop over nodes of the element

    sum_Kp.inv(inv_sum_Kp);

    MD.m_res_at_point.fill(0.0);

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      MD.m_res_at_point += MD.m_grad_F[dim].row_transpose(q);
    }

    for (Uint m = 0; m < MD.nb_nodes(); ++m)
    {
      // I) Accumulate into residuals (local RHS)
      math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(m * NEQ, NEQ);
      //      vector += scalar * (   matrix   *  matrix   * vector )
      //      [NEQ]  += scalar * ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      // res_in_node += wj_q * (MD.m_Kp[m] * MD.inv_sum_Kp *
      // MD.m_res_at_point);

      beta_distr_mat = Kp[m] * inv_sum_Kp;
      res_in_node += wj_q * (beta_distr_mat * MD.m_res_at_point);

      // II) Accumulate into local Jacobi matrix (local LHS)
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
template <typename MD1, typename MD2>
void PGLDAImplementation<Physics>::compute_adv_residuals_no_flux_metric(
    math::DenseConstMatView<Real> const &sol_nodal_values,
    interpolation::CellGeoMetric<MD1> const &CGM, interpolation::CellSolutionMetric<MD2> const &CSM,
    typename Physics::FluxV const &cell_residual, RDMethodScratchData<Physics> &MD)
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename CellGeoMetType, typename CellSolMetType, typename CellFluxMetType>
void PGLDAImplementation<Physics>::compute_adv_reaction_residuals(
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

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA
  math::DenseSVec<Real, NEQ> eigvalues;
#endif

  /// ASSEMBLE THE NODAL RESIDUALS
  using JM = typename phys_model::JM;
  /*
  JM &Rv = MD.Rv();
  JM &Lv = MD.Lv();
  JM &Dvp = MD.Dvp();
  JM &Dvm = MD.Dvm();
  */
  JM &sum_Kp     = MD.sum_Kp();
  JM &inv_sum_Kp = MD.inv_sum_Kp();

  // Vector of matrices K(+), one matrix per element DOF
  common::ArrayView<typename phys_model::JM, _1D, Uint> Kp = MD.Kp();

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
        MD.m_grad_sf_at_pt_ref[dim] = MD.m_dV_u[dim](q, n);
      }

      MD.m_grad_sf_at_pt_phys = inv_J * MD.m_grad_sf_at_pt_ref;

#if USE_SIMPLIFIED_EVALUATION_OF_K_PGLDA
      Physics::build_K_mat(MD.m_props, MD.m_grad_sf_at_pt_phys, eigvalues, Kp[n],
                           [](const Real eigvalue) { return (std::max(0.0, eigvalue)); });

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < NEQ; ++e)
      {
        MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, std::abs(eigvalues[e]));
      }

      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

#else

      Physics::flux_jacobian_eigen_structure(MD.m_props, MD.m_grad_sf_at_pt_phys, Rv, Lv, Dvp);

      MD.m_max_eigenvalue = 1.e-6;

      for (Uint e = 0; e < Physics::NEQ; ++e)
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

      // MD.m_elem_wave_speed[n] += wj_q * (MD.m_max_eigenvalue);
      MD.m_elem_wave_speed[n] = std::max(MD.m_elem_wave_speed[n], MD.m_max_eigenvalue);

      m_Kp[n] = Rv * Dvp * Lv;
#endif

      sum_Kp += Kp[n];
    }

    sum_Kp.inv(inv_sum_Kp);

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

    for (Uint n = 0; n < MD.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MD.elem_node_res_block(n * NEQ, NEQ);
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      res_in_node += wj_q * (Kp[n] * inv_sum_Kp * MD.m_res_at_point);
      // MD.elem_node_res[n] += wj_q * (MD.m_Kp[n] * MD.inv_sum_Kp *
      // MD.m_res_at_point) -
      //                          wj_q * MD.m_V(q, n) *
      //                          sq.row_transpose(q);
    }

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
struct CellSchemeSelector<Physics, PGLDA>
{
  typedef PGLDAImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
