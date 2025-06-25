#ifndef PDEKIT_RDM_Solver_FACET_LDA_hpp
#define PDEKIT_RDM_Solver_FACET_LDA_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

template <typename Physics>
class FacetLDA
{
  public:
  enum
  {
    needs_volume_derivatives_on_trace = 1
  };

  /// TYPEDEFS
  typedef Physics phys_model;

  /// Constructor
  FacetLDA();

  /// Destructor
  ~FacetLDA();

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, Uint ScratchDataDim>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_left,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_right,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
                                 mesh::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
                                 mesh::QuadraturePermutation const quad_perm_right,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDRight);

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD3 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
                                 interpolation::CellGeoMetric<MD2> const &CellLeftGM,
                                 interpolation::CellGeoMetric<MD2> const &CellRightGM,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_left,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_right,
                                 interpolation::CellSolutionMetric<MD3> const &CellLeftContourSM,
                                 mesh::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD3> const &CellRightContourSM,
                                 mesh::QuadraturePermutation const quad_perm_right,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDRight);

  private:
};

// ----------------------------------------------------------------------------

template <typename Physics>
FacetLDA<Physics>::FacetLDA()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
FacetLDA<Physics>::~FacetLDA()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, Uint ScratchDataDim>
void FacetLDA<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
    math::DenseConstMatView<Real> const &facet_nodal_values_left,
    math::DenseConstMatView<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
    mesh::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
    mesh::QuadraturePermutation const quad_perm_right,
    RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
    RDMethodScratchData<Physics, ScratchDataDim> &MDRight)

{
  // Reset the element residuals
  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n].fill(0.0);
  }

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n].fill(0.0);
  }

  // Reset the element update coefficient
  MDRight.m_elem_wave_speed.fill(0.0);
  MDLeft.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &w = FacetGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = FacetGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq = FacetGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Facet normals at quadrature points
  const math::DenseConstMatView<Real> facet_normals = FacetGM.normals();
  // std::cout << "Facet normals " << std::endl << facet_normals << std::endl;

  // Solution in quadrature points
  const math::DenseConstMatView<Real> uq_left = FacetSMLeft.field_values();
  // std::cout << "Solution at quadrature points (left) = " << std::endl <<
  // uq_left << std::endl;
  const math::DenseConstMatView<Real> uq_right = FacetSMRight.field_values();
  // std::cout << "Solution at quadrature points (right) = " << std::endl <<
  // uq_right << std::endl;

  /// ASSEMBLE THE NODAL RESIDUALS
  Real elem_volume = 0.0;

  for (Uint q = 0; q < FacetGM.nb_qd_pts(); ++q)
  {
    const Uint qL = q; // quad_perm_left.get().vertex(q);
    const Uint qR = q; // quad_perm_right.get().vertex(q);

    MDLeft.sum_Kp.fill(0.0);
    MDRight.sum_Kp.fill(0.0);

    const Real wj_q = w[qL] * jdet[qL];
    elem_volume += wj_q;

    const math::DenseConstVecView<Real> normal_q = facet_normals.row_transpose(qL);

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      MDLeft.m_grad_u_at_point.insert_row(eq, normal_q);
      MDRight.m_grad_u_at_point.insert_row(eq, normal_q);
    }

    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      for (Uint d = 0; d < Physics::DIM; ++d)
      {
        MDLeft.m_grad_u_at_point(eq, d) *= 1.0;
      }
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(q), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(Xq.row_transpose(q), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    Physics::flux_jacobian_eigen_structure(MDLeft.m_props, normal_q, MDLeft.Rv, MDLeft.Lv,
                                           MDLeft.Dvp);

    Physics::flux_jacobian_eigen_structure(MDRight.m_props, -1. * normal_q, MDRight.Rv, MDRight.Lv,
                                           MDRight.Dvp);

    Physics::flux(MDLeft.m_props, normal_q, MDLeft.m_flux_integral_in_elem);
    Physics::flux(MDRight.m_props, normal_q, MDRight.m_flux_integral_in_elem);

    MDLeft.m_res_at_point.fill(0.0);
    MDLeft.m_res_at_point = MDRight.m_flux_integral_in_elem - MDLeft.m_flux_integral_in_elem;

    // const math::ConstMatrixBlock<Real> inv_J = FacetGM.inv_jacobi(q);

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      MDLeft.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, std::abs(MDLeft.Dvp(e, e)));
        MDLeft.Dvp(e, e)        = std::max(1.e-8, MDLeft.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      MDLeft.m_elem_wave_speed[n] += wj_q * MDLeft.m_max_eigenvalue;

      MDLeft.m_Kp[n] = MDLeft.Rv * MDLeft.Dvp * MDLeft.Lv;
      // Here we accumulate to the sum(K+) of the LEFT data
      MDLeft.sum_Kp += MDLeft.m_Kp[n];
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      MDRight.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(MDRight.Dvp(e, e)));
        MDRight.Dvp(e, e)        = std::max(1.e-8, MDRight.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      MDRight.m_elem_wave_speed[n] += wj_q * MDRight.m_max_eigenvalue;

      MDRight.m_Kp[n] = MDRight.Rv * MDRight.Dvp * MDRight.Lv;
      // Here we accumulate to the sum(K+) of the LEFT data AGAIN!
      MDLeft.sum_Kp += MDRight.m_Kp[n];
    }

    MDLeft.sum_Kp.inv(MDLeft.inv_sum_Kp);

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      MDLeft.m_elem_node_res[n] +=
          wj_q * (MDLeft.m_Kp[n] * MDLeft.inv_sum_Kp * MDLeft.m_res_at_point);
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      MDRight.m_elem_node_res[n] +=
          wj_q * (MDRight.m_Kp[n] * MDLeft.inv_sum_Kp * MDLeft.m_res_at_point);
    }

  } // Loop over quadrature points

  /*
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_update_coeff[n] /= elem_volume;
  }
  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_update_coeff[n] /= elem_volume;
  }
  */
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
void FacetLDA<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
    interpolation::CellGeoMetric<MD2> const &CellLeftGM,
    interpolation::CellGeoMetric<MD2> const &CellRightGM,
    math::DenseConstMatView<Real> const &facet_nodal_values_left,
    math::DenseConstMatView<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD3> const &CellLeftContourSM,
    mesh::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD3> const &CellRightContourSM,
    mesh::QuadraturePermutation const quad_perm_right,
    RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
    RDMethodScratchData<Physics, ScratchDataDim> &MDRight)
{
  // Reset the element residuals
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n].fill(0.0);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n].fill(0.0);
  }

  // Reset the element update coefficient
  MDLeft.m_elem_wave_speed.fill(0.0);
  MDRight.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &w = FacetGM.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Quadrature weights - left cell
  const math::DenseDVec<Real> &wL = CellLeftGM.pt_weights();

  // Quadrature weights - right cell
  const math::DenseDVec<Real> &wR = CellRightGM.pt_weights();

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdet = FacetGM.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Jacobian determinants - left cell
  const math::DenseConstVecView<Real> jdetL = CellLeftGM.jdet();

  // Jacobian determinants - right cell
  const math::DenseConstVecView<Real> jdetR = CellRightGM.jdet();

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> Xq_facet = FacetGM.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  const math::DenseConstMatView<Real> Xq_cell_left  = CellLeftGM.interpolated_coords();
  const math::DenseConstMatView<Real> Xq_cell_right = CellRightGM.interpolated_coords();

  // Facet normals at quadrature points
  const math::DenseConstMatView<Real> facet_normals = FacetGM.normals();
  // std::cout << "Facet normals " << std::endl << facet_normals << std::endl;

  // Solution in quadrature points
  const math::DenseConstMatView<Real> uq_left = CellLeftContourSM.field_values();
  // std::cout << "Solution at quadrature points (left) = " << std::endl <<
  // uq_left << std::endl;
  const math::DenseConstMatView<Real> uq_right = CellRightContourSM.field_values();
  // std::cout << "Solution at quadrature points (right) = " << std::endl <<
  // uq_right << std::endl;

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MDLeft.m_grad_u[d]  = CellLeftContourSM.field_derivatives(d);
    MDRight.m_grad_u[d] = CellRightContourSM.field_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  /*
  std::cout << "Left subset = ";
  for (Uint i = 0; i < MDLeft.m_sub_idx.size(); ++i)
  {
    std::cout << MDLeft.m_sub_idx[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Right subset = ";
  for (Uint i = 0; i < MDRight.m_sub_idx.size(); ++i)
  {
    std::cout << MDRight.m_sub_idx[i] << " ";
  }
  std::cout << std::endl << std::endl;
  */

  /*
  std::cout << "Nodes left = " << MDLeft.m_nb_nodes << std::endl;
  std::cout << "Nodes right = " << MDRight.m_nb_nodes << std::endl;
  */

  /*
    std::cout << "Quad left = " << quad_perm_left.get().code() << std::endl;
    std::cout << "Quad right = " << quad_perm_right.get().code() << std::endl;
    for (Uint q = 0; q < FacetGM.nb_qd_pts(); ++q)
    {
      const Uint qL = quad_perm_left.get().vertex(q);
      const Uint qR = quad_perm_right.get().vertex(q);
      std::cout << "[" << qL << "," << qR << "] = [" << Xq_cell_left.row(qL)
    << "  "
                << Xq_cell_right.row(qR) << "]" << std::endl;
      for(Uint c = 0; c < Xq_cell_left.cols(); ++c)
      {
        if ( std::abs(Xq_cell_left(qL,c) - Xq_cell_right(qR,c)) > 1.e-10 )
        {
          std::cout << "Values differ [qL,qR] = [" << qL << "," << qR << "]"
    << std::endl;
        }
      }
    }
    std::cout << "**************************" << std::endl;
    */

  /*
  for (Uint q = 0; q < FacetGM.nb_qd_pts(); ++q)
  {
    const Uint qL = quad_perm_left.get().vertex(q);
    const Uint qR = quad_perm_right.get().vertex(q);
    std::cout << "[" << qL << "," << qR << "] = [" << uq_left.row(qL) << "  "
  << uq_right.row(qR)
              << "]" << std::endl;
  }
  std::cout << "**************************" << std::endl;
  */

  // std::cout << "Coordinates on facet:" << std::endl;
  // std::cout << Xq_facet << std::endl;
  // std::cout << "Coordinates left adjacent cell:" << std::endl;
  // std::cout << Xq_cell_left << std::endl;
  // std::cout << "Coordinates right adjacent cell:" << std::endl;
  // std::cout << Xq_cell_right << std::endl;
  // std::cout << "*****************************************************" <<
  // std::endl;

  /// ASSEMBLE THE NODAL RESIDUALS
  for (Uint q = 0; q < FacetGM.nb_qd_pts(); ++q)
  {
    MDLeft.sum_Kp.fill(0.0);
    MDRight.sum_Kp.fill(0.0);

    const Uint qL = quad_perm_left.get().vertex(q);
    const Uint qR = quad_perm_right.get().vertex(q);

    const Real wj_qL = w[q] * jdet[q]; // wL[qL] * jdetL[qL];
    const Real wj_qR = w[q] * jdet[q]; // wR[qR] * jdetR[qR];

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MDLeft.m_grad_u_at_point.insert_col(dim, MDLeft.m_grad_u[dim].row(qL));
      MDRight.m_grad_u_at_point.insert_col(dim, MDRight.m_grad_u[dim].row(qR));
    }

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq_cell_left.row_transpose(qL), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(Xq_cell_right.row_transpose(qR), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    const math::DenseConstMatView<Real> inv_J_L = CellLeftGM.inv_jacobi(qL);
    const math::DenseConstMatView<Real> inv_J_R = CellRightGM.inv_jacobi(qR);

    const math::DenseConstVecView<Real> normal_q = facet_normals.row_transpose(q); // SHOULD BE q?

    Physics::flux(MDLeft.m_props, normal_q, MDLeft.m_flux_integral_in_elem);
    Physics::flux(MDRight.m_props, normal_q, MDRight.m_flux_integral_in_elem);

    math::DenseSVec<Real, Physics::NEQ> eig_left, eig_right;
    eig_left.fill(0.0);
    eig_right.fill(0.0);

    Physics::flux_jacobian_eigen_values(MDLeft.m_props, normal_q, eig_left);
    Physics::flux_jacobian_eigen_values(MDRight.m_props, -1. * normal_q, eig_right);

    Real flux_max_eig_value = 1.e-8;

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MDLeft.m_grad_sf_at_pt_ref[dim] = MDLeft.m_dV_u[dim](qL, n);
      }

      MDLeft.m_grad_sf_at_pt_phys = inv_J_L * MDLeft.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MDLeft.m_props, MDLeft.m_grad_sf_at_pt_phys, MDLeft.Rv,
                                             MDLeft.Lv, MDLeft.Dvp);

      MDLeft.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, std::abs(MDLeft.Dvp(e, e)));
        flux_max_eig_value      = std::max(flux_max_eig_value, std::abs(eig_left[e]));
        MDLeft.Dvp(e, e)        = std::max(1.e-8, MDLeft.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      // MDLeft.m_elem_update_coeff[n] += wj_qL * MDLeft.m_max_eigenvalue;
      MDLeft.m_elem_wave_speed[n] += wL[qL] * jdetL[qL] * MDLeft.m_max_eigenvalue;

      MDLeft.m_Kp[n] = MDLeft.Rv * MDLeft.Dvp * MDLeft.Lv;
      MDLeft.sum_Kp += MDLeft.m_Kp[n];
    }

    MDLeft.sum_Kp.inv(MDLeft.inv_sum_Kp);

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MDRight.m_grad_sf_at_pt_ref[dim] = MDRight.m_dV_u[dim](qR, n);
      }

      MDRight.m_grad_sf_at_pt_phys = inv_J_R * MDRight.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MDRight.m_props, MDRight.m_grad_sf_at_pt_phys,
                                             MDRight.Rv, MDRight.Lv, MDRight.Dvp);

      MDRight.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < Physics::NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the
        // Jacobian???) Option a: spectral radius
        MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(MDRight.Dvp(e, e)));
        flux_max_eig_value       = std::max(flux_max_eig_value, std::abs(eig_right[e]));
        MDRight.Dvp(e, e)        = std::max(1.e-8, MDRight.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e,
        // e));
      }

      // MDRight.m_elem_update_coeff[n] += wj_qR *
      // MDRight.m_max_eigenvalue;
      MDRight.m_elem_wave_speed[n] += wR[qR] * jdetR[qR] * MDRight.m_max_eigenvalue;

      MDRight.m_Kp[n] = MDRight.Rv * MDRight.Dvp * MDRight.Lv;
      MDRight.sum_Kp += MDRight.m_Kp[n];
    }

    MDRight.sum_Kp.inv(MDRight.inv_sum_Kp);

    // This is numerical flux at given quadrature point
    // MDRight.m_res_at_point.fill(0.0);
    MDRight.m_res_at_point =
        0.5 * (MDLeft.m_flux_integral_in_elem + MDRight.m_flux_integral_in_elem) -
        0.5 * flux_max_eig_value * (uq_right.row_transpose(qR) - uq_left.row_transpose(qL));

#if 0
    math::StaticVector<Real, Physics::NEQ> tmp;

    if ((-1.1 <= Xq_cell_left(0, X0)) && (Xq_cell_left(0, X0) <= -0.9) && (Xq_cell_left(0, X1) <= 0.5))
    {
      std::cout << std::endl;
      std::cout << "---------------------------" << std::endl;
      std::cout << "f* @ q[" << q << "] = " << MDRight.m_res_at_point << std::endl << std::endl;
      std::cout << "Eig left = " << eig_left << std::endl;
      std::cout << "Eig right = " << eig_right << std::endl;
      std::cout << "Eig max = " << flux_max_eig_value << std::endl;
      std::cout << "---------------------------" << std::endl;
      std::cout << std::endl;
      std::cout << "fL @ qL[" << q << "] = " << MDLeft.m_flux_integral_in_elem << std::endl;
      tmp = MDRight.m_res_at_point - MDLeft.m_flux_integral_in_elem;
      std::cout << "(f*-fL) @ qL[" << q << "] = " << tmp << std::endl;
      for(Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
      {
        std::cout << "K+[" << n << "] @ qL[" << q << "] = " << MDLeft.m_Kp[n];
      }

      std::cout << std::endl;
      std::cout << "fR @ qR[" << q << "] = " << MDRight.m_flux_integral_in_elem << std::endl;
      tmp = MDRight.m_res_at_point - MDRight.m_flux_integral_in_elem;
      std::cout << "(f*-fR) @ qR[" << q << "] = " << tmp << std::endl;
      for(Uint n = 0; n < MDRight.m_nb_nodes; ++n)
      {
        std::cout << "K+[" << n << "] @ qR[" << q << "] = " << MDRight.m_Kp[n];
      }
    }
#endif

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      MDLeft.m_elem_node_res[n] +=
          wj_qL * (MDLeft.m_Kp[n] * MDLeft.inv_sum_Kp *
                   (MDRight.m_res_at_point - MDLeft.m_flux_integral_in_elem));
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      MDRight.m_elem_node_res[n] -=
          wj_qR * (MDRight.m_Kp[n] * MDRight.inv_sum_Kp *
                   (MDRight.m_res_at_point - MDRight.m_flux_integral_in_elem));
    }

  } // Loop over quadrature points

#if 0
  if ((-1.1 <= Xq_cell_left(0, X0)) && (Xq_cell_left(0, X0) <= -0.9) && (Xq_cell_left(0, X1) <= 0.5))
  {
    std::cout << "Coords left = " << std::endl << Xq_cell_left;
    std::cout << "(local idx =";
    for(Uint i = 0; i < MDLeft.m_sub_idx.size(); ++i)
    {
      std::cout << " " << MDLeft.m_sub_idx[i];
    }
    std::cout << ")" << std::endl << std::endl;
    std::cout << "Coords right = " << std::endl << Xq_cell_right;
    std::cout << "(local idx =";
    for(Uint i = 0; i < MDRight.m_sub_idx.size(); ++i)
    {
      std::cout << " " << MDRight.m_sub_idx[i];
    }
    std::cout << ")" << std::endl << std::endl;

    std::cout << "uqL = " << std::endl;
    std::cout << uq_left << std::endl;

    std::cout << "uqR = " << std::endl;
    std::cout << uq_right << std::endl;

    std::cout << "Residuals left = " << std::endl;
    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      std::cout << MDLeft.m_elem_node_res[n] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Residuals right = " << std::endl;
    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      std::cout << MDRight.m_elem_node_res[n] << std::endl;
    }

    std::cout << "***********************************************************" << std::endl;
  }
#endif
}

// ----------------------------------------------------------------------------

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
