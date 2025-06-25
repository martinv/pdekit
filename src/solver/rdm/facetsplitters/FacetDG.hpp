#ifndef PDEKIT_Solver_RDM_FACET_DG_hpp
#define PDEKIT_Solver_RDM_FACET_DG_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"
#include "solver/NumFlux.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/facetsplitters/FacetSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Tag class to mark LDA with variable beta coefficients

// ----------------------------------------------------------------------------

class FacetDG
{
};

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics>
class FacetDGImplementation
{
  public:
  enum
  {
    needs_volume_derivatives_on_trace = 0
  };

  /// TYPEDEFS
  typedef Physics phys_model;
  typedef RDMethodScratchData<Physics, Physics::DIM - 1> method_data;

  /// Constructor
  FacetDGImplementation();

  /// Destructor
  ~FacetDGImplementation();

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1> const &FacetGMLeft,
                                 interpolation::CellGeoMetricWithNormals<MD1> const &FacetGMRight,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_left,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_right,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
                                 mesh::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
                                 mesh::QuadraturePermutation const quad_perm_right,
                                 interpolation::FacetFluxMetric<MD3, Physics> const &FacetFMLeft,
                                 interpolation::FacetFluxMetric<MD3, Physics> const &FacetFMRight,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
                                 RDMethodScratchData<Physics, ScratchDataDim> &MDRight);

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD3 ... type of metric data for the cell metric terms

  /*
  template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1>
  const &FacetGMLeft, interpolation::CellGeoMetricWithNormals<MD2> const
  &FacetGMRight, interpolation::CellGeoMetric<MD2> const &CellLeftGM,
                                 interpolation::CellGeoMetric<MD2> const
  &CellRightGM, math::ConstMatrixBlock<Real> const &facet_nodal_values_left,
                                 math::ConstMatrixBlock<Real> const
  &facet_nodal_values_right, interpolation::CellSolutionMetric<MD3> const
  &FacetSMLeft, quadrature::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD3> const
  &FacetSMRight, quadrature::QuadraturePermutation const quad_perm_right,
                                 RDMethodScratchData<Physics, ScratchDataDim>
  &MDLeft, RDMethodScratchData<Physics, ScratchDataDim> &MDRight);
  */

  private:
  enum
  {
    NEQ = phys_model::NEQ
  };

  NumFluxLaxFriedrichs<Physics> m_num_flux;
  //  NumFluxAUSM<Physics> m_num_flux;
};

// ----------------------------------------------------------------------------

template <typename Physics>
FacetDGImplementation<Physics>::FacetDGImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
FacetDGImplementation<Physics>::~FacetDGImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
void FacetDGImplementation<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGMLeft,
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGMRight,
    math::DenseConstMatView<Real> const &facet_nodal_values_left,
    math::DenseConstMatView<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
    mesh::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
    mesh::QuadraturePermutation const quad_perm_right,
    interpolation::FacetFluxMetric<MD3, Physics> const &FacetFMLeft,
    interpolation::FacetFluxMetric<MD3, Physics> const &FacetFMRight,
    RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
    RDMethodScratchData<Physics, ScratchDataDim> &MDRight)

{
  // Reset the element residuals
  MDLeft.elem_node_res().fill(0.0);
  MDRight.elem_node_res().fill(0.0);

  // Reset the element update coefficient
  MDLeft.m_elem_wave_speed.fill(0.0);
  MDRight.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DenseDVec<Real> &wL = FacetGMLeft.pt_weights();
  const math::DenseDVec<Real> &wR = FacetGMRight.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::DenseConstVecView<Real> jdetL = FacetGMLeft.jdet();
  const math::DenseConstVecView<Real> jdetR = FacetGMRight.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet <<
  // std::endl;

  // Coordinates of integration points in physical space
  const math::DenseConstMatView<Real> XqL = FacetGMLeft.interpolated_coords();
  const math::DenseConstMatView<Real> XqR = FacetGMRight.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Facet normals at quadrature points
  const math::DenseConstMatView<Real> facet_normalsL = FacetGMLeft.normals();
  const math::DenseConstMatView<Real> facet_normalsR = FacetGMRight.normals();
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

  MDLeft.m_flux_integral_in_elem.fill(0.0);
  MDRight.m_flux_integral_in_elem.fill(0.0);

  // std::cout << "***********************************" << std::endl;

  /*
  std::cout << "Nq/geo = " << FacetGM.nb_qd_pts() << std::endl;
  std::cout << "Nq/sol = " << quad_perm_left.get().size() << std::endl;
  */

  // std::cout << "************************************************" <<
  // std::endl;
  for (Uint qG = 0; qG < FacetGMLeft.nb_qd_pts(); ++qG)
  {
    /*
    const Uint qL = q; // quad_perm_left.get().vertex(q);
    const Uint qR = q; // quad_perm_right.get().vertex(q);
    */

    const Uint qL = quad_perm_left.get().vertex(qG);
    const Uint qR = quad_perm_right.get().vertex(qG);

    // std::cout << "qL = " << qL << ", qR = " << qR << std::endl;

    const Real wj_qL = wL[qL] * jdetL[qL];
    const Real wj_qR = wR[qR] * jdetR[qR];
    elem_volume += wj_qL;

    const math::DenseConstVecView<Real> normal_qL = facet_normalsL.row_transpose(qL);
    const math::DenseConstVecView<Real> normal_qR = facet_normalsR.row_transpose(qR);

#if 0
    /*
    std::cout << "nL = " << normal_qL << std::endl;
    std::cout << "nR = " << normal_qR << std::endl;
    std::cout << "xL = " << XqL.row_transpose(qL) << std::endl;
    std::cout << "xR = " << XqR.row_transpose(qR) << std::endl << std::endl;
    */

    /// Compute properties at each quadrature point
    Physics::compute_properties(XqL.row_transpose(qL), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(XqR.row_transpose(qR), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    Physics::flux(MDLeft.m_props, normal_qL, MDLeft.m_norm_flux_at_point);
    Physics::flux(MDRight.m_props, normal_qR, MDRight.m_norm_flux_at_point);

    MDLeft.m_flux_integral_in_elem += wj_qL * MDLeft.m_norm_flux_at_point;
    MDRight.m_flux_integral_in_elem += wj_qR * MDRight.m_norm_flux_at_point;

    math::StaticVector<Real, NEQ> eig_left, eig_right;
    eig_left.fill(0.0);
    eig_right.fill(0.0);

    Physics::flux_jacobian_eigen_values(MDLeft.m_props, normal_qL, eig_left);
    Physics::flux_jacobian_eigen_values(MDRight.m_props, normal_qR, eig_right);

    MDLeft.m_max_eigenvalue = 1.e-8;

    for (Uint e = 0; e < NEQ; ++e)
    {
      // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
      // __positive__ eigenvalues, or should it be the biggest
      // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
      // Option a: spectral radius
      MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, std::abs(eig_left[e]));
      eig_left[e] = std::max(1.e-8, eig_left[e]); // THIS IS REDUNDANT

      // Option b: largest positive eigenvalue
      // MDLeft.Dvp(e, e) = std::max(0.0, MDLeft.Dvp(e, e));
      // MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, MDLeft.Dvp(e, e));
    }

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      // MDLeft.m_elem_update_coeff[n] += wj_q * MDLeft.m_max_eigenvalue;
      MDLeft.m_elem_wave_speed[n] = std::max(MDLeft.m_elem_wave_speed[n], MDLeft.m_max_eigenvalue);
    }

    MDRight.m_max_eigenvalue = 1.e-8;

    for (Uint e = 0; e < NEQ; ++e)
    {
      // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
      // __positive__ eigenvalues, or should it be the biggest
      // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
      // Option a: spectral radius
      MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(eig_right[e]));
      eig_right[e] = std::max(1.e-8, eig_right[e]); // THIS IS REDUNDANT

      // Option b: largest positive eigenvalue
      // MDRight.Dvp(e, e) = std::max(1.e-8, MDRight.Dvp(e, e));
      // MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, MDRight.Dvp(e, e));
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      // MDRight.m_elem_update_coeff[n] += wj_q * MDRight.m_max_eigenvalue;
      MDRight.m_elem_wave_speed[n] =
          std::max(MDRight.m_elem_wave_speed[n], MDRight.m_max_eigenvalue);
    }

    // This is numerical flux at given quadrature point
    MDRight.m_res_at_point = 0.5 * (MDLeft.m_norm_flux_at_point - MDRight.m_norm_flux_at_point) -
                             0.5 * std::max(MDLeft.m_max_eigenvalue, MDRight.m_max_eigenvalue) *
                                 (uq_right.row_transpose(qR) - uq_left.row_transpose(qL));

#else

    m_num_flux.compute(XqL.row_transpose(qL), XqR.row_transpose(qR), normal_qL,
                       uq_left.row_transpose(qL), uq_right.row_transpose(qR),
                       MDRight.m_res_at_point);

    const Real max_eig_LR =
        std::max(m_num_flux.max_eigvalue_left(), m_num_flux.max_eigvalue_right());

    for (Uint n = 0; n < MDLeft.nb_nodes(); ++n)
    {
      // MDLeft.m_elem_wave_speed[n] =
      //     std::max(MDLeft.m_elem_wave_speed[n],
      //     m_num_flux.max_eigvalue_left());
      MDLeft.m_elem_wave_speed[n] = std::max(MDLeft.m_elem_wave_speed[n], max_eig_LR);
    }
    for (Uint n = 0; n < MDRight.nb_nodes(); ++n)
    {
      // MDRight.m_elem_wave_speed[n] =
      //    std::max(MDRight.m_elem_wave_speed[n],
      //    m_num_flux.max_eigvalue_right());
      MDRight.m_elem_wave_speed[n] = std::max(MDRight.m_elem_wave_speed[n], max_eig_LR);
    }

    /*
    MDLeft.m_norm_flux_at_point = m_num_flux.flux_left();
    MDRight.m_norm_flux_at_point = -1. * m_num_flux.flux_right();
    */

    // Compute the normal flux coming from the left facet
    MDLeft.m_norm_flux_at_point.fill(0.0);
    for (Uint d = 0; d < phys_model::DIM; ++d)
    {
      MDLeft.m_norm_flux_at_point += normal_qL[d] * FacetFMLeft.flux_values(d).row_transpose(qL);
    }

    // Compute the normal flux coming from the right facet
    MDRight.m_norm_flux_at_point.fill(0.0);
    for (Uint d = 0; d < phys_model::DIM; ++d)
    {
      MDRight.m_norm_flux_at_point += normal_qR[d] * FacetFMRight.flux_values(d).row_transpose(qR);
    }

#endif

    for (Uint n = 0; n < MDLeft.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MDLeft.elem_node_res_block(n * NEQ, NEQ);
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      res_in_node +=
          wj_qL * MDLeft.m_V_u(qL, n) * (MDRight.m_res_at_point - MDLeft.m_norm_flux_at_point);
    }

    for (Uint n = 0; n < MDRight.nb_nodes(); ++n)
    {
      math::DenseVecView<Real> res_in_node = MDRight.elem_node_res_block(n * NEQ, NEQ);
      //            vector += scalar *  (   matrix   *  matrix   * vector
      //            ) [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]
      //            )
      res_in_node -=
          wj_qR * MDRight.m_V_u(qR, n) * (MDRight.m_res_at_point + MDRight.m_norm_flux_at_point);
    }

  } // Loop over quadrature points

  const Real wave_speed_factor =
      ((FacetSMLeft.std_region_type().cell_transform_id() == mesh::CellTransform::NO_TRANS) &&
       FacetSMRight.std_region_type().cell_transform_id() == mesh::CellTransform::NO_TRANS)
          ? 1.0
          : 2.0;

  for (Uint n = 0; n < MDLeft.nb_nodes(); ++n)
  {
    MDLeft.m_elem_wave_speed[n] *= wave_speed_factor;
  }
  for (Uint n = 0; n < MDRight.nb_nodes(); ++n)
  {
    MDRight.m_elem_wave_speed[n] *= wave_speed_factor;
    ;
  }

  /*
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    // MDLeft.m_elem_wave_speed[n] *= elem_volume;
    MDLeft.m_elem_wave_speed[n] *= elem_volume * elem_volume;
  }
  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    // MDRight.m_elem_wave_speed[n] *= elem_volume;
    MDRight.m_elem_wave_speed[n] *= elem_volume * elem_volume;
  }
  */
}

// ----------------------------------------------------------------------------

#if 0
template <typename Physics>
template <typename MD1, typename MD2, typename MD3, Uint ScratchDataDim>
void FacetDGImplementation<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGMLeft,
    interpolation::CellGeoMetricWithNormals<MD2> const &FacetGMRight,
    interpolation::CellGeoMetric<MD2> const &CellLeftGM,
    interpolation::CellGeoMetric<MD2> const &CellRightGM,
    math::ConstMatrixBlock<Real> const &facet_nodal_values_left,
    math::ConstMatrixBlock<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD3> const &FacetSMLeft,
    quadrature::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD3> const &FacetSMRight,
    quadrature::QuadraturePermutation const quad_perm_right,
    RDMethodScratchData<Physics, ScratchDataDim> &MDLeft,
    RDMethodScratchData<Physics, ScratchDataDim> &MDRight)
{
  // Reset the element residuals
  MDRight.elem_node_res().fill(0.0);
  MDLeft.elem_node_res().fill(0.0);

  // Reset the element update coefficient
  MDRight.m_elem_wave_speed.fill(0.0);
  MDLeft.m_elem_wave_speed.fill(0.0);

  // Quadrature weights
  const math::DynamicVector<Real> &w = FacetGMLeft.pt_weights();
  // std::cout << "Quadrature weights:" << std::endl << w << std::endl;

  // Jacobian determinants in quadrature points
  const math::ConstVectorBlock<Real> jdet = FacetGMLeft.jdet();
  // std::cout << "Jacobian determinants = " << std::endl << jdet << std::endl;

  // Coordinates of integration points in physical space
  const math::ConstMatrixBlock<Real> Xq = FacetGMLeft.interpolated_coords();
  // std::cout << "Quadrature points = " << std::endl << Xq << std::endl;

  // Facet normals at quadrature points
  const math::ConstMatrixBlock<Real> facet_normals = FacetGMLeft.normals();
  // std::cout << "Facet normals " << std::endl << facet_normals << std::endl;

  // Solution in quadrature points
  const math::ConstMatrixBlock<Real> uq_left = FacetSMLeft.field_values();
  // std::cout << "Solution at quadrature points (left) = " << std::endl << uq_left << std::endl;
  const math::ConstMatrixBlock<Real> uq_right = FacetSMRight.field_values();
  // std::cout << "Solution at quadrature points (right) = " << std::endl << uq_right << std::endl;

  // Gradients of solution
  // std::cout << "Solution gradients at quadrature points = " << std::endl;
  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    // std::cout << "d = " << d << std::endl;
    MDRight.m_grad_u[d] = FacetSMRight.field_derivatives(d);
    MDLeft.m_grad_u[d] = FacetSMLeft.field_derivatives(d);
    // std::cout << m_grad_u[d] << std::endl;
  }

  /// ASSEMBLE THE NODAL RESIDUALS
  Real elem_volume = 0.0;

  MDLeft.m_flux_integral_in_elem.fill(0.0);
  MDRight.m_flux_integral_in_elem.fill(0.0);

  for (Uint qG = 0; qG < FacetGMLeft.nb_qd_pts(); ++qG)
  {
    const Uint qL = quad_perm_left.get().vertex(qG);
    const Uint qR = quad_perm_right.get().vertex(qG);

    MDRight.sum_Kp.fill(0.0);
    MDLeft.sum_Kp.fill(0.0);

    const Real wj_q = w[qL] * jdet[qL];
    elem_volume += wj_q;

    for (Uint dim = 0; dim < Physics::DIM; ++dim)
    {
      // Gradient of solution
      // The size of MD.m_grad_u_at_point = [NEQ x DIM]
      MDLeft.m_grad_u_at_point.insert_col(dim, MDLeft.m_grad_u[dim].row(qL));
      MDRight.m_grad_u_at_point.insert_col(dim, MDRight.m_grad_u[dim].row(qR));
    }

    /// Compute properties at each quadrature point
    // Note that the physical coordinate Xq is taken at qL! It has to coincide
    // with the physical coordinate seen from the 'left' side of the facet ...
    Physics::compute_properties(Xq.row_transpose(qG), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(Xq.row_transpose(qG), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    const math::ConstMatrixBlock<Real> inv_J_L = CellLeftGM.inv_jacobi(qL);
    const math::ConstMatrixBlock<Real> inv_J_R = CellRightGM.inv_jacobi(qR);

    const math::ConstVectorBlock<Real> normal_q = facet_normals.row_transpose(qL);

    Physics::flux(MDLeft.m_props, normal_q, MDLeft.m_norm_flux_at_point);
    Physics::flux(MDRight.m_props, normal_q, MDRight.m_norm_flux_at_point);

    MDLeft.m_flux_integral_in_elem += wj_q * MDLeft.m_norm_flux_at_point;
    MDRight.m_flux_integral_in_elem += wj_q * MDRight.m_norm_flux_at_point;

    MDRight.m_res_at_point.fill(0.0);
    MDRight.m_res_at_point = MDRight.m_flux_integral_in_elem - MDLeft.m_flux_integral_in_elem;


    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MDRight.m_grad_sf_at_pt_ref[dim] = MDRight.m_dV[dim](qR, n);
      }

      MDRight.m_grad_sf_at_pt_phys = inv_J_R * MDRight.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MDRight.m_props, MDRight.m_grad_sf_at_pt_phys,
                                             MDRight.Rv, MDRight.Lv, MDRight.Dvp);

      MDRight.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
        // Option a: spectral radius
        MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(MDRight.Dvp(e, e)));
        MDRight.Dvp(e, e) = std::max(1.e-8, MDRight.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e, e));
      }

      // MDRight.m_elem_update_coeff[n] += wj_q * (MDRight.m_max_eigenvalue);
      MDRight.m_elem_wave_speed[n] +=
          std::max(MDRight.m_elem_wave_speed[n], MDRight.m_max_eigenvalue);

      MDRight.m_Kp[n] = MDRight.Rv * MDRight.Dvp * MDRight.Lv;
      MDRight.sum_Kp += MDRight.m_Kp[n];
    }

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      /// Gradient of solution shape functions
      for (Uint dim = 0; dim < Physics::DIM; ++dim)
      {
        // [DIM][dim] = (nb_qd_pts,nb_nodes)(q,n)
        MDLeft.m_grad_sf_at_pt_ref[dim] = MDLeft.m_dV[dim](qL, n);
      }

      MDLeft.m_grad_sf_at_pt_phys = inv_J_L * MDLeft.m_grad_sf_at_pt_ref;

      Physics::flux_jacobian_eigen_structure(MDLeft.m_props, MDLeft.m_grad_sf_at_pt_phys, MDLeft.Rv,
                                             MDLeft.Lv, MDLeft.Dvp);

      MDLeft.m_max_eigenvalue = 1.e-8;

      for (Uint e = 0; e < NEQ; ++e)
      {
        // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
        // __positive__ eigenvalues, or should it be the biggest
        // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
        // Option a: spectral radius
        MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, std::abs(MDLeft.Dvp(e, e)));
        MDLeft.Dvp(e, e) = std::max(1.e-8, MDLeft.Dvp(e, e));

        // Option b: largest positive eigenvalue
        // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));
        // MD.m_max_eigenvalue = std::max(MD.m_max_eigenvalue, MD.Dvp(e, e));
      }

      // MDLeft.m_elem_update_coeff[n] += wj_q * (MDLeft.m_max_eigenvalue);
      MDLeft.m_elem_wave_speed[n] += std::max(MDLeft.m_elem_wave_speed[n], MDLeft.m_max_eigenvalue);

      MDLeft.m_Kp[n] = MDLeft.Rv * MDLeft.Dvp * MDLeft.Lv;
      // Here we accumulate to the sum(K+) of the RIGHT data
      MDRight.sum_Kp += MDLeft.m_Kp[n];
    }

    MDRight.sum_Kp.inv(MDRight.inv_sum_Kp);

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      math::DenseVecView<Real> res_in_node = MDRight.elem_node_res_block(n * NEQ, NEQ);

      //   vector += scalar *  (   matrix   *  matrix   * vector )
      //   [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (MDRight.m_Kp[n] * MDRight.inv_sum_Kp * MDRight.m_res_at_point);
    }

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      math::DenseVecView<Real> res_in_node = MDLeft.elem_node_res_block(n * NEQ, NEQ);

      //   vector += scalar *  (   matrix   *  matrix   * vector )
      //   [NEQ]  += scalar *  ( [NEQxNEQ]  * [NEQxNEQ] * [NEQ]  )
      res_in_node += wj_q * (MDLeft.m_Kp[n] * MDRight.inv_sum_Kp * MDRight.m_res_at_point);
    }

  } // Loop over quadrature points


  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_wave_speed[n] = 0.0; //*= elem_volume;
  }
  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_wave_speed[n] = 0.0; //*= elem_volume;
  }
}
#endif

// ----------------------------------------------------------------------------

template <typename Physics>
struct FacetSchemeSelector<Physics, FacetDG>
{
  typedef FacetDGImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
