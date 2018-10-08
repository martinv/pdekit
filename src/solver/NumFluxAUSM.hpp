#ifndef PDEKIT_Solver_Num_Flux_AUSM_hpp
#define PDEKIT_Solver_Num_Flux_AUSM_hpp

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
class NumFluxAUSM
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
  /// METHODS

  static Real M1_plus(const Real M);
  static Real M1_minus(const Real M);

  static Real M2_plus(const Real M);
  static Real M2_minus(const Real M);

  static Real p3_plus(const Real M);
  static Real p3_minus(const Real M);

  /// DATA

  enum
  {
    NEQ = Physics::NEQ
  };

  typename Physics::Properties m_props_L;
  typename Physics::Properties m_props_R;
  typename Physics::FluxV m_flux_L;
  typename Physics::FluxV m_flux_R;
  typename Physics::FluxV m_conv_flux;

  typename Physics::Properties::SolGradM m_grad_u;

  math::DenseSVec<Real, NEQ> m_eig_left;
  math::DenseSVec<Real, NEQ> m_eig_right;

  Real m_eig_max_L;
  Real m_eig_max_R;
};

// ----------------------------------------------------------------------------

template <typename Physics>
template <typename V1, typename V2, typename V3, typename V4, typename V5, typename V6>
void NumFluxAUSM<Physics>::compute(const V1 &x_L, const V2 &x_R, const V3 &normal, const V4 &u_L,
                                   const V5 &u_R, V6 &flux)
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

  Real vn_L = 0.0;
  Real vn_R = 0.0;

  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    vn_L += m_props_L.V[d] * normal[d];
    vn_R += m_props_R.V[d] * normal[d];
  }

  const Real M_L = vn_L / m_props_L.a;
  const Real M_R = vn_R / m_props_R.a;

  const Real M_half = M2_plus(M_L) + M2_minus(M_R);

  // Convective part
  if (M_half >= 0.0)
  {
    for (Uint eq = 0; eq < (NEQ - 1); ++eq)
    {
      flux[eq] = M_half * m_props_L.a * m_props_L.rho;
    }
    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      flux[d + 1] *= m_props_L.V[d];
    }

    flux[NEQ - 1] = M_half * m_props_L.a * (m_props_L.rhoE + m_props_L.P);
  }
  else
  {
    for (Uint eq = 0; eq < (NEQ - 1); ++eq)
    {
      flux[eq] = M_half * m_props_R.a * m_props_R.rho;
    }
    for (Uint d = 0; d < Physics::DIM; ++d)
    {
      flux[d + 1] *= m_props_R.V[d];
    }
    flux[NEQ - 1] = M_half * m_props_R.a * (m_props_R.rhoE + m_props_R.P);
  }

  // Add pressure part
  const Real p_half = p3_plus(M_L) * m_props_L.P + p3_minus(M_R) * m_props_R.P;

  for (Uint d = 0; d < Physics::DIM; ++d)
  {
    flux[d + 1] += p_half * normal[d];
  }
}

// ----------------------------------------------------------------------------

template <typename Physics>
Real NumFluxAUSM<Physics>::max_eigvalue_left() const
{
  return m_eig_max_L;
}

// ----------------------------------------------------------------------------

template <typename Physics>
Real NumFluxAUSM<Physics>::max_eigvalue_right() const
{
  return m_eig_max_R;
}

// ----------------------------------------------------------------------------

template <typename Physics>
const typename Physics::FluxV &NumFluxAUSM<Physics>::flux_left() const
{
  return m_flux_L;
}

// ----------------------------------------------------------------------------

template <typename Physics>
const typename Physics::FluxV &NumFluxAUSM<Physics>::flux_right() const
{
  return m_flux_R;
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::M1_plus(const Real M)
{
  return 0.5 * (M + std::abs(M));
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::M1_minus(const Real M)
{
  return 0.5 * (M - std::abs(M));
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::M2_plus(const Real M)
{
  if (std::abs(M) <= 1.)
  {
    return 0.25 * (M + 1.) * (M + 1.);
  }
  return 0.5 * (M + std::abs(M)); // This is M1_plus(M)
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::M2_minus(const Real M)
{
  if (std::abs(M) <= 1.)
  {
    return -0.25 * (M - 1.) * (M - 1.);
  }
  return 0.5 * (M - std::abs(M)); // This is M1_minus(M)
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::p3_plus(const Real M)
{
  const Real abs_M = std::abs(M);
  if (abs_M <= 1.)
  {
    return 0.25 * (M + 1.) * (M + 1.) * (2. - M); // M2_plus(M) * (2-M)
  }
  // Supersonic case:
  return 0.5 / M * (M + abs_M);
}

// ----------------------------------------------------------------------------

template <typename Physics>
inline Real NumFluxAUSM<Physics>::p3_minus(const Real M)
{
  const Real abs_M = std::abs(M);
  if (abs_M <= 1.)
  {
    return 0.25 * (M - 1.) * (M - 1.) * (2. + M); // - M2_minus(M) * (2+M)
  }
  // Supersonic case:
  return 0.5 / M * (M - abs_M);
}

// ----------------------------------------------------------------------------

} // namespace solver

} // namespace pdekit

#endif
