#ifndef PDEKIT_RDM_Solver_FACET_LF_hpp
#define PDEKIT_RDM_Solver_FACET_LF_hpp

#include "interpolation/CellFluxMetric.hpp"
#include "interpolation/CellGeoMetric.hpp"
#include "interpolation/CellSolutionMetric.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"
#include "solver/rdm/RDMethodScratchData.hpp"
#include "solver/rdm/facetsplitters/FacetSchemeSelector.hpp"

namespace pdekit
{

namespace solver
{

namespace rdm
{

// Tag class to mark Lax-Friedrichs scheme for facets

// ----------------------------------------------------------------------------

class FacetLF
{
};

// ----------------------------------------------------------------------------

namespace internal
{

template <typename Physics, Uint NEQ>
struct LFResidualLimiter;

template <typename Physics>
class FacetLFImplementation
{
  public:
  enum
  {
    needs_volume_derivatives_on_trace = 0
  };

  /// TYPEDEFS
  typedef Physics phys_model;
  typedef PGRDLFMethodData<Physics, Physics::DIM - 1> method_data;

  /// Constructor
  FacetLFImplementation();

  /// Destructor
  ~FacetLFImplementation();

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD2 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_left,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_right,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
                                 mesh::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
                                 mesh::QuadraturePermutation const quad_perm_right,
                                 RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
                                 RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight);

  /// Compute the residuals on one element
  /// @param MD1 ... type of metric data for the facet metric terms
  /// @param MD3 ... type of metric data for the cell metric terms
  template <typename MD1, typename MD2, typename MD3>
  void facet_nodal_contributions(interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
                                 interpolation::CellGeoMetric<MD2> const &CellLeftGM,
                                 interpolation::CellGeoMetric<MD2> const &CellRightGM,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_left,
                                 math::DenseConstMatView<Real> const &facet_nodal_values_right,
                                 interpolation::CellSolutionMetric<MD3> const &FacetSMLeft,
                                 mesh::QuadraturePermutation const quad_perm_left,
                                 interpolation::CellSolutionMetric<MD3> const &FacetSMRight,
                                 mesh::QuadraturePermutation const quad_perm_right,
                                 RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
                                 RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight);

  private:
  /// TYPES
  typedef LFResidualLimiter<Physics, Physics::NEQ> residual_limiter_type;

  /// Maximum eigenvalues of the decomposed jacobian \grad(F) \cdot \nabla
  /// \varphi where F is the tensor of inviscid fluxes
  Real m_alpha;

  /// Normalized velocity vector
  typename Physics::FluxV m_norm_velocity;

  /// Average state
  math::DenseSVec<Real, Physics::NEQ> m_u_avg;

  /// Integral of the jump in normal fluxes
  typename Physics::FluxV m_facet_residual;

  /// Sum of x^K+ or x^{\Gamma}+ contributions for all dofs in the element
  typename Physics::FluxV m_sum_xp;
};

// ----------------------------------------------------------------------------

template <typename Physics>
FacetLFImplementation<Physics>::FacetLFImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
FacetLFImplementation<Physics>::~FacetLFImplementation()
{
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2>
void FacetLFImplementation<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
    math::DenseConstMatView<Real> const &facet_nodal_values_left,
    math::DenseConstMatView<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD2> const &FacetSMLeft,
    mesh::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD2> const &FacetSMRight,
    mesh::QuadraturePermutation const quad_perm_right,
    RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
    RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight)

{
  // Reset the facet residuals
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n].fill(0.0);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n].fill(0.0);
  }

  // Reset the facet update coefficient
  MDLeft.m_elem_wave_speed.fill(0.0);
  MDRight.m_elem_wave_speed.fill(0.0);

  m_u_avg.fill(0.0);

  // facet_nodal_values_{left|right} is a matrix of dimension
  // [(nb. nodes) x (nb. equations)]

  const Uint tot_nb_states_on_facet          = MDLeft.m_nb_nodes + MDRight.m_nb_nodes;
  const Real one_over_tot_nb_states_on_facet = 1. / tot_nb_states_on_facet;

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    m_u_avg += facet_nodal_values_left.row_transpose(n);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    m_u_avg += facet_nodal_values_right.row_transpose(n);
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    m_u_avg[eq] *= one_over_tot_nb_states_on_facet;
  }

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

  m_facet_residual.fill(0.0);

  m_alpha = 1.e-9;

  // Real dot_prod = 0.0;

  MDLeft.m_flux_integral_in_elem.fill(0.0);
  MDRight.m_flux_integral_in_elem.fill(0.0);

  for (Uint qG = 0; qG < FacetGM.nb_qd_pts(); ++qG)
  {
    /*
    const Uint qL = q; // quad_perm_left.get().vertex(q);
    const Uint qR = q; // quad_perm_right.get().vertex(q);
    */

    const Uint qL = quad_perm_left.get().vertex(qG);
    const Uint qR = quad_perm_right.get().vertex(qG);

    // std::cout << "qL = " << qL << ", qR = " << qR << std::endl;

    const Real wj_q = w[qL] * jdet[qL];

    const math::DenseConstVecView<Real> normal_q = facet_normals.row_transpose(qL);

    /// Compute properties at each quadrature point
    Physics::compute_properties(Xq.row_transpose(qG), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(Xq.row_transpose(qG), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    Physics::flux(MDLeft.m_props, normal_q, MDLeft.m_norm_flux_at_point);
    Physics::flux(MDRight.m_props, normal_q, MDRight.m_norm_flux_at_point);

    // Compute the gradient of each shape function in physical space
    // const math::ConstMatrixBlock<Real> inv_J_L =
    // CellLeftGM.inv_jacobi(qL); const math::ConstMatrixBlock<Real> inv_J_R
    // = CellRightGM.inv_jacobi(qR);

    // m_facet_flux_jump_integral += wj_q * (MDRight.m_norm_flux_at_point -
    // MDLeft.m_norm_flux_at_point);

    MDLeft.m_flux_integral_in_elem += wj_q * MDLeft.m_norm_flux_at_point;
    MDRight.m_flux_integral_in_elem += wj_q * MDRight.m_norm_flux_at_point;

    Physics::flux_jacobian_eigen_structure(MDLeft.m_props, MDLeft.m_props.V, MDLeft.Rv, MDLeft.Lv,
                                           MDLeft.Dvp);

    Physics::flux_jacobian_eigen_structure(MDRight.m_props, MDRight.m_props.V, MDRight.Rv,
                                           MDRight.Lv, MDRight.Dvp);

    MDLeft.m_max_eigenvalue = 1.e-8;

    for (Uint e = 0; e < Physics::NEQ; ++e)
    {
      // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
      // __positive__ eigenvalues, or should it be the biggest
      // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
      // Option a: spectral radius
      MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue, std::abs(MDLeft.Dvp(e, e)));
    }

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      // MDLeft.m_elem_wave_speed[n] += wj_q * MDLeft.m_max_eigenvalue;
      MDLeft.m_elem_wave_speed[n] = std::max(MDLeft.m_elem_wave_speed[n], MDLeft.m_max_eigenvalue);
    }

    MDRight.m_max_eigenvalue = 1.e-8;

    for (Uint e = 0; e < Physics::NEQ; ++e)
    {
      // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
      // __positive__ eigenvalues, or should it be the biggest
      // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
      // Option a: spectral radius
      MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(MDRight.Dvp(e, e)));
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      // MDRight.m_elem_wave_speed[n] += wj_q * MDRight.m_max_eigenvalue;
      MDRight.m_elem_wave_speed[n] =
          std::max(MDRight.m_elem_wave_speed[n], MDRight.m_max_eigenvalue);
    }

  } // Loop over quadrature points

  // -------------------------------------
  // Integral of flux jump along the facet
  // -------------------------------------
  // The facet residual is in fact an integral of the jump of the fluxes along
  // the facet res(facet) = integral ( f_R - f_L )
  m_facet_residual = MDRight.m_flux_integral_in_elem - MDLeft.m_flux_integral_in_elem;

  // -------------
  // Compute alpha
  // -------------

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    m_alpha = std::max(std::abs(MDLeft.m_elem_wave_speed[n]), m_alpha);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    m_alpha = std::max(std::abs(MDRight.m_elem_wave_speed[n]), m_alpha);
  }

  // ---------------------------
  // Compute nodal residuals
  // ---------------------------

  MDLeft.m_elem_node_res[0]  = 1.0 / tot_nb_states_on_facet * m_facet_residual;
  MDRight.m_elem_node_res[0] = MDLeft.m_elem_node_res[0];

  for (Uint n = 1; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n]  = MDLeft.m_elem_node_res[0];
    MDRight.m_elem_node_res[n] = MDLeft.m_elem_node_res[0];
  }

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n] += m_alpha * (facet_nodal_values_left.row_transpose(n) - m_u_avg);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n] += m_alpha * (facet_nodal_values_right.row_transpose(n) - m_u_avg);
  }

  // Limiting for HO scheme
  residual_limiter_type::limit_facet_residuals(m_facet_residual, m_sum_xp, MDLeft, MDRight);
}

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename MD1, typename MD2, typename MD3>
void FacetLFImplementation<Physics>::facet_nodal_contributions(
    interpolation::CellGeoMetricWithNormals<MD1> const &FacetGM,
    interpolation::CellGeoMetric<MD2> const &CellLeftGM,
    interpolation::CellGeoMetric<MD2> const &CellRightGM,
    math::DenseConstMatView<Real> const &facet_nodal_values_left,
    math::DenseConstMatView<Real> const &facet_nodal_values_right,
    interpolation::CellSolutionMetric<MD3> const &FacetSMLeft,
    mesh::QuadraturePermutation const quad_perm_left,
    interpolation::CellSolutionMetric<MD3> const &FacetSMRight,
    mesh::QuadraturePermutation const quad_perm_right,
    RDMethodScratchData<Physics, Physics::DIM - 1> &MDLeft,
    RDMethodScratchData<Physics, Physics::DIM - 1> &MDRight)
{
  // Reset the facet residuals
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n].fill(0.0);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n].fill(0.0);
  }

  // Reset the facet update coefficient
  MDLeft.m_elem_wave_speed.fill(0.0);
  MDRight.m_elem_wave_speed.fill(0.0);

  m_u_avg.fill(0.0);

  // facet_nodal_values_{left|right} is a matrix of dimension
  // [(nb. nodes) x (nb. equations)]

  const Uint tot_nb_states_on_facet          = MDLeft.m_nb_nodes + MDRight.m_nb_nodes;
  const Real one_over_tot_nb_states_on_facet = 1. / tot_nb_states_on_facet;

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    m_u_avg += facet_nodal_values_left.row_transpose(n);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    m_u_avg += facet_nodal_values_right.row_transpose(n);
  }

  for (Uint eq = 0; eq < Physics::NEQ; ++eq)
  {
    m_u_avg[eq] *= one_over_tot_nb_states_on_facet;
  }

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

  m_facet_residual.fill(0.0);

  m_alpha = 1.e-9;

  // Real dot_prod = 0.0;

  MDLeft.m_flux_integral_in_elem.fill(0.0);
  MDRight.m_flux_integral_in_elem.fill(0.0);

  for (Uint qG = 0; qG < FacetGM.nb_qd_pts(); ++qG)
  {
    const Uint qL = quad_perm_left.get().vertex(qG);
    const Uint qR = quad_perm_right.get().vertex(qG);

    Physics::compute_properties(Xq.row_transpose(qG), uq_left.row_transpose(qL),
                                MDLeft.m_grad_u_at_point, MDLeft.m_props);

    Physics::compute_properties(Xq.row_transpose(qG), uq_right.row_transpose(qR),
                                MDRight.m_grad_u_at_point, MDRight.m_props);

    const math::DenseConstVecView<Real> normal = facet_normals.row_transpose(qL);

    Physics::flux(MDLeft.m_props, normal, MDLeft.m_norm_flux_at_point);
    Physics::flux(MDRight.m_props, normal, MDRight.m_norm_flux_at_point);

    const Real wj_q = w[qL] * jdet[qL];

    MDLeft.m_flux_integral_in_elem += wj_q * MDLeft.m_norm_flux_at_point;
    MDRight.m_flux_integral_in_elem += wj_q * MDRight.m_norm_flux_at_point;

    Physics::flux_jacobian_eigen_structure(MDLeft.m_props, MDLeft.m_props.V, MDLeft.Rv, MDLeft.Lv,
                                           MDLeft.Dvp);

    Physics::flux_jacobian_eigen_structure(MDRight.m_props, MDRight.m_props.V, MDRight.Rv,
                                           MDRight.Lv, MDRight.Dvp);

    /*
    for (Uint eq = 0; eq < Physics::NEQ; ++eq)
    {
      const Real eig_plus_left = std::max(std::abs(MDLeft.Dvp(eq, eq)),
    m_alpha); const Real eig_plus_right = std::max(std::abs(MDRight.Dvp(eq,
    eq)), m_alpha); m_alpha = std::max(eig_plus_left, eig_plus_right);

      const Real vn_left = MDLeft.m_props.V[X] * normal[X] +
    MDLeft.m_props.V[Y]
    * normal[Y]; const Real vn_right = MDRight.m_props.V[X] * normal[X] +
    MDRight.m_props.V[Y] * normal[Y];

      const Real eig_plus_left = std::max(std::abs(vn_left), m_alpha);
      const Real eig_plus_right = std::max(std::abs(vn_right), m_alpha);

      m_alpha = std::max(m_alpha, std::abs(vn_left - vn_right));
    }
    */

    MDLeft.m_max_eigenvalue = MDRight.m_max_eigenvalue = 1.e-6;

    for (Uint e = 0; e < Physics::NEQ; ++e)
    {
      // FIXME: SHOULD m_max_eigenvalue be simply the biggest of all
      // __positive__ eigenvalues, or should it be the biggest
      // eigenvalue in magnitude (i.e. spectral radius of the Jacobian???)
      // Option a: spectral radius
      MDLeft.m_max_eigenvalue  = std::max(MDLeft.m_max_eigenvalue, std::abs(MDLeft.Dvp(e, e)));
      MDRight.m_max_eigenvalue = std::max(MDRight.m_max_eigenvalue, std::abs(MDRight.Dvp(e, e)));
      // MD.Dvp(e, e) = std::max(0.0, MD.Dvp(e, e));

      // THIS IS NOT CORRECT
      // Option b: largest positive eigenvalue
      // MDLeft.Dvp(e, e) = std::max(0.0, MDLeft.Dvp(e, e));
      // MDRight.Dvp(e, e) = std::max(0.0, MDRight.Dvp(e, e));
      // MDLeft.m_max_eigenvalue = std::max(MDLeft.m_max_eigenvalue,
      // MDLeft.Dvp(e, e)); MDRight.m_max_eigenvalue =
      // std::max(MDRight.m_max_eigenvalue, MDRight.Dvp(e, e));
    }

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      MDLeft.m_elem_wave_speed[n] +=
          one_over_tot_nb_states_on_facet * wj_q * MDLeft.m_max_eigenvalue;
    }
    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      MDRight.m_elem_wave_speed[n] +=
          one_over_tot_nb_states_on_facet * wj_q * MDRight.m_max_eigenvalue;
    }

    /*
    const Real vn_left = MDLeft.m_props.V[X] * normal[X] +
    MDLeft.m_props.V[Y] * normal[Y]; const Real vn_right =
    MDRight.m_props.V[X] * normal[X] + MDRight.m_props.V[Y] * normal[Y];

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      // MDLeft.m_elem_update_coeff[n] += one_over_tot_nb_states_on_facet *
    wj_q
    * std::abs(vn_right
      // - vn_left);
      MDLeft.m_elem_wave_speed[n] +=
          one_over_tot_nb_states_on_facet * wj_q *
    std::max(std::abs(vn_left), std::abs(vn_right));
    }

    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      // MDRight.m_elem_update_coeff[n] += one_over_tot_nb_states_on_facet *
    wj_q *
      // std::abs(vn_right - vn_left);
      MDRight.m_elem_wave_speed[n] +=
          one_over_tot_nb_states_on_facet * wj_q *
    std::max(std::abs(vn_left), std::abs(vn_right));
    }
    */

  } // Loop over quadrature points

  // -------------------------------------
  // Integral of flux jump along the facet
  // -------------------------------------
  // The facet residual is in fact an integral of the jump of the fluxes along
  // the facet res(facet) = integral ( f_R - f_L )
  m_facet_residual = MDRight.m_flux_integral_in_elem - MDLeft.m_flux_integral_in_elem;

  // -------------
  // Compute alpha
  // -------------

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    m_alpha = std::max(std::abs(MDLeft.m_elem_wave_speed[n]), m_alpha);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    m_alpha = std::max(std::abs(MDRight.m_elem_wave_speed[n]), m_alpha);
  }

  // ---------------------------
  // Compute nodal residuals
  // ---------------------------

  MDLeft.m_elem_node_res[0]  = 1.0 / tot_nb_states_on_facet * m_facet_residual;
  MDRight.m_elem_node_res[0] = MDLeft.m_elem_node_res[0];

  for (Uint n = 1; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n]  = MDLeft.m_elem_node_res[0];
    MDRight.m_elem_node_res[n] = MDLeft.m_elem_node_res[0];
  }

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_node_res[n] += m_alpha * (facet_nodal_values_left.row_transpose(n) - m_u_avg);
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_node_res[n] += m_alpha * (facet_nodal_values_right.row_transpose(n) - m_u_avg);
  }

  // ---------------------------
  // Compute update coefficients
  // ---------------------------

  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_wave_speed[n] = m_alpha - MDLeft.m_elem_wave_speed[n]; // / MD.m_nb_nodes;
  }

  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_wave_speed[n] = m_alpha - MDRight.m_elem_wave_speed[n]; // / MD.m_nb_nodes;
  }

  /*
  for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
  {
    MDLeft.m_elem_update_coeff[n] = m_alpha / tot_nb_states_on_facet;
  }
  for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
  {
    MDRight.m_elem_update_coeff[n] = m_alpha / tot_nb_states_on_facet;
  }
  */

  /*
  for (Uint q = 0; q < CGM.nb_qd_pts(); ++q)
  {
    const Real wj_q = w[q] * jdet[q];

    for (Uint n = 0; n < MDLeft.m_nb_nodes; ++n)
    {
      MDLeft.m_elem_update_coeff[n] += wj_q * (m_max_eigenvalue);
    }
    for (Uint n = 0; n < MDRight.m_nb_nodes; ++n)
    {
      MDRight.m_elem_update_coeff[n] += wj_q * (m_max_eigenvalue);
    }
  }
  */

  // Limiting for HO scheme
  residual_limiter_type::limit_facet_residuals(m_facet_residual, m_sum_xp, MDLeft, MDRight);
}

// ----------------------------------------------------------------------------

template <typename Physics>
struct FacetSchemeSelector<Physics, FacetLF>
{
  typedef FacetLFImplementation<Physics> type;
};

// ----------------------------------------------------------------------------

} // namespace internal

} // namespace rdm

} // namespace solver

} // namespace pdekit

#endif
