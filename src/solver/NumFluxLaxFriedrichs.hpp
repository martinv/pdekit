#ifndef PDEKIT_Solver_Num_Flux_LaxFriedrichs_hpp
#define PDEKIT_Solver_Num_Flux_LaxFriedrichs_hpp

#include "common/Meta.hpp"
#include "common/PDEKit.hpp"
#include "math/DenseSVec.hpp"
#include "math/TensorRank.hpp"

namespace pdekit
{

namespace solver
{

// ----------------------------------------------------------------------------

template <typename Physics>
class NumFluxLaxFriedrichs
{
  public:
  template <typename V1, typename V2, typename V3, typename V4, typename V5, typename V6>
  void compute(const V1 &x_L, const V2 &x_R, const V3 &normal, const V4 &u_L, const V5 &u_R,
               V6 &flux);

  Real max_eigvalue_left() const;

  Real max_eigvalue_right() const;

  const typename Physics::FluxV &flux_left() const;

  const typename Physics::FluxV &flux_right() const;

  private:
  enum
  {
    NEQ = Physics::NEQ
  };

  typename Physics::Properties m_props_L;
  typename Physics::Properties m_props_R;
  typename Physics::FluxV m_flux_L;
  typename Physics::FluxV m_flux_R;
  typename Physics::Properties::SolGradM m_grad_u;

  math::DenseSVec<Real, NEQ> m_eig_left;
  math::DenseSVec<Real, NEQ> m_eig_right;

  Real m_eig_max_L;
  Real m_eig_max_R;
};

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename V1, typename V2, typename V3, typename V4, typename V5, typename V6>
void NumFluxLaxFriedrichs<Physics>::compute(const V1 &x_L, const V2 &x_R, const V3 &normal,
                                            const V4 &u_L, const V5 &u_R, V6 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'x_L' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'x_R' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V3>::value, math::tensor_rank_1>::value,
                "'normal' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V4>::value, math::tensor_rank_1>::value,
                "'u_L' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V5>::value, math::tensor_rank_1>::value,
                "'u_R' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V6>::value, math::tensor_rank_1>::value,
                "'flux' must be a vector");

  Physics::compute_properties(x_L, u_L, m_grad_u, m_props_L);
  Physics::compute_properties(x_R, u_R, m_grad_u, m_props_R);

  Physics::flux(m_props_L, normal, m_flux_L);
  Physics::flux(m_props_R, normal, m_flux_R);

  m_eig_left.fill(0.0);
  m_eig_right.fill(0.0);

  Physics::flux_jacobian_eigen_values(m_props_L, normal, m_eig_left);
  Physics::flux_jacobian_eigen_values(m_props_R, normal, m_eig_right);

  m_eig_max_L = 1.e-8;
  m_eig_max_R = 1.e-8;

  for (Uint e = 0; e < NEQ; ++e)
  {
    // m_eig_left is the biggest eigenvalue in magnitude (i.e. spectral
    // radius of the Jacobian
    m_eig_max_L = std::max(m_eig_max_L, std::abs(m_eig_left[e]));
    // m_eig_left[e] = std::max(1.e-8, m_eig_right[e]); // This is redundant

    m_eig_max_R = std::max(m_eig_max_R, std::abs(m_eig_right[e]));
    // m_eig_right[e] = std::max(1.e-8, m_eig_right[e]); // This is
    // redundant
  }

  // This is the resulting numerical flux
  flux = 0.5 * (m_flux_L + m_flux_R) - 0.5 * std::max(m_eig_max_L, m_eig_max_R) * (u_R - u_L);
}

// ----------------------------------------------------------------------------

template <typename Physics>
Real NumFluxLaxFriedrichs<Physics>::max_eigvalue_left() const
{
  return m_eig_max_L;
}

// ----------------------------------------------------------------------------

template <typename Physics>
Real NumFluxLaxFriedrichs<Physics>::max_eigvalue_right() const
{
  return m_eig_max_R;
}

// ----------------------------------------------------------------------------

template <typename Physics>
const typename Physics::FluxV &NumFluxLaxFriedrichs<Physics>::flux_left() const
{
  return m_flux_L;
}

// ----------------------------------------------------------------------------

template <typename Physics>
const typename Physics::FluxV &NumFluxLaxFriedrichs<Physics>::flux_right() const
{
  return m_flux_R;
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
