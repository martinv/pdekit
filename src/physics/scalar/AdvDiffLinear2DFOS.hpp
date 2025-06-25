#ifndef PDEKIT_Physics_Adv_Diff_Linear2D_FOS_hpp
#define PDEKIT_Physics_Adv_Diff_Linear2D_FOS_hpp

#include <array>

#include "physics/PhysModelT.hpp"
#include "physics/scalar/AdvDiff2DFOSProperties.hpp"

namespace pdekit
{

namespace physics
{

class AdvDiffLinear2DFOS : public PhysModelT<AdvDiff2DFOSProperties>
{

  public:
  /// TYPEDEFS

  using Properties    = typename PhysModelT<AdvDiff2DFOSProperties>::Properties;
  using physics_class = PhysicsClassScalar;

  /// Constructor
  AdvDiffLinear2DFOS();

  /// Destructor
  ~AdvDiffLinear2DFOS();

  // --------------------------------------------------------------------------

  /// Compute the physical properties at given point
  /// @param coord     ... point coordinates in which the properties are
  /// computed
  /// @param sol       ... solution vector 'U' (as input for properties
  /// computation)
  /// @param grad_vars ... gradient of 'U'
  /// @return p        ... computed properties

  template <typename V1, typename V2, typename M1>
  static void compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars,
                                 Properties &p);

  // --------------------------------------------------------------------------

  /// Compute the value of the solution 'u' from the properties
  /// @param p     ... properties
  /// @return vars ... computed variables

  template <typename V1>
  static void compute_variables(const Properties &p, V1 &vars);

  // --------------------------------------------------------------------------

  /// Compute the fluxes
  /// @param p     ... input physical properties to compute the flux
  /// @return flux ... matrix of fluxes computed from the properties
  template <typename M1>
  static void flux(const Properties &p, M1 &flux);

  // --------------------------------------------------------------------------

  /// Compute the fluxes in given direction
  /// @param p         ... input properties to compute the flux
  /// @param direction ... direction vector in which the fluxes should be
  /// computed
  /// @return flux     ... resulting flux, flux = Flux_x * direction[X] +
  /// Flux_y
  /// * direction[Y]
  template <typename V1, typename V2>
  static void flux(const Properties &p, const V1 &direction, V2 &flux);

  // --------------------------------------------------------------------------

  /// Compute the source term
  /// @param p          ... input properties to compute the flux
  /// @return src_term  ... computed source term
  template <typename V1>
  static void source_term(const Properties &p, V1 &src_term);

  // --------------------------------------------------------------------------

  /// Compute the eigenvalues of the flux jacobian
  /// @param p         ... physical properties
  /// @param direction ... direction vector, the jacobian that we
  ///                      decompose will be J = direction[X]*dFx/dU +
  /// direction[Y]*dFy/dU
  /// @return Dv       ... vector of computed eigenvalues
  template <typename V1, typename V2>
  static void flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv);

  // --------------------------------------------------------------------------

  /// Compute the eigen values of the flux jacobians
  /// @param p         ... physical properties
  /// @param direction ... direction vector, the Jacobian that we
  ///                      decompose will be J = direction[X]*dFx/dU +
  /// direction[Y]*dFy/dU
  /// @return Dv       ... vector of computed eigenvalues
  /// @param  Op       ... extra operator applied to the eigen values (max,
  /// '+' operator for example)
  template <typename V1, typename V2, typename OP>
  static void flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv, OP &op);

  // --------------------------------------------------------------------------

  /// Decompose the eigen structure of the flux jacobians projected on the
  /// gradients
  /// @param p         ... physical properties
  /// @param direction ... direction onto which the Jacobian
  ///                      is projected (J = direction[X]*dFx/dU +
  /// direction[Y]*dFy/dU)
  /// @return Rv       ... matrix of right eigenvectors
  /// @return Lv       ... matrix of left eigenvectors
  /// @return Dv       ... matrix of eigenvalues
  template <typename V1, typename M1, typename M2, typename M3>
  static void flux_jacobian_eigen_structure(const Properties &p, const V1 &direction, M1 &Rv,
                                            M2 &Lv, M3 &Dv);
  // --------------------------------------------------------------------------

  /// Decompose the eigen structure of the flux jacobians projected on the
  /// gradients
  /// @param p           ... physical properties
  /// @param direction   ... direction onto which the Jacobian
  ///                               is projected (J = direction[X]*dFx/dU +
  ///                               direction[Y]*dFy/dU)
  /// @return eigvalues  ... vector of eigenvalues of the matrix K^{+}
  /// @return K_mat      ... matrix K^{+}

  template <typename V1, typename V2, typename M1, typename OP>
  static void build_K_mat(const Properties &p, const V1 &direction, V2 &eigvalues, M1 &K_mat,
                          const OP &op);

  // --------------------------------------------------------------------------

  /// compute the PDE residual
  /// @param p          ... physical properties
  /// @param flux_jacob ... vector of flux jacobians
  /// @return res       ... computed residual (lambda_x * Jx + lambda_y * Jy)
  ///                       (the advection vector is lambda = (lambda_x,
  /// lambda_y)
  template <typename M1, typename V1>
  static void residual(const Properties &p, std::array<M1, DIM> &flux_jacob, V1 &res);
};

// --------------------------------------------------------------------------

template <typename V1, typename V2, typename M1>
void AdvDiffLinear2DFOS::compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars,
                                            Properties &p)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'coord' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'sol' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'grad_vars' must be a matrix");

  p.coords    = coord;     // cache the coordinates locally
  p.vars      = sol;       // cache the variables locally
  p.grad_vars = grad_vars; // cache the gradient of variables locally

  p.V[X] = 1; // Linear advection
  p.V[Y] = 0; //

  p.mu = 1.0; // Constant diffusion
  p.Lr = 0.5;
}

// ----------------------------------------------------------------------------

template <typename V1>
void AdvDiffLinear2DFOS::compute_variables(const Properties &p, V1 &vars)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'vars' must be a vector");
  vars = p.vars;
}

// ----------------------------------------------------------------------------

template <typename M1>
void AdvDiffLinear2DFOS::flux(const Properties &p, M1 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'flux' must be a matrix");
  const Real inv_Tr = p.mu / (p.Lr * p.Lr);

  flux(0, X0) = -p.mu * p.vars[1];
  flux(1, X0) = -p.vars[0] * inv_Tr;
  flux(2, X0) = 0.0;

  flux(0, X1) = -p.mu * p.vars[2];
  flux(1, X1) = 0.0;
  flux(2, X1) = -p.vars[0] * inv_Tr;
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2>
void AdvDiffLinear2DFOS::flux(const Properties &p, const V1 &direction, V2 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'flux' must be a vector");

  const Real inv_Tr = p.mu / (p.Lr * p.Lr);
  flux[0]           = -p.mu * (direction[X0] * p.vars[1] + direction[X1] * p.vars[2]);
  flux[1]           = -direction[0] * p.vars[0] * inv_Tr;
  flux[2]           = -direction[1] * p.vars[0] * inv_Tr;
}

// ----------------------------------------------------------------------------

template <typename V1>
inline void AdvDiffLinear2DFOS::source_term(const Properties &p, V1 &src_term)
{
  const Real inv_Tr = p.mu / (p.Lr * p.Lr);
  src_term[0]       = 0.0;
  src_term[1]       = p.vars[1] * inv_Tr;
  src_term[2]       = p.vars[2] * inv_Tr;
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2>
void AdvDiffLinear2DFOS::flux_jacobian_eigen_values(const Properties &p, const V1 &direction,
                                                    V2 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");

  const Real inv_Tr = p.mu / (p.Lr * p.Lr);
  const Real norm2  = direction[0] * direction[0] + direction[1] * direction[1];

  Dv[0] = -std::sqrt(norm2 * p.mu * inv_Tr);
  Dv[1] = Dv[0];
  Dv[2] = 0.0;
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2, typename OP>
void AdvDiffLinear2DFOS::flux_jacobian_eigen_values(const Properties &p, const V1 &direction,
                                                    V2 &Dv, OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");

  const Real inv_Tr     = p.mu / (p.Lr * p.Lr);
  const Real norm2      = direction[0] * direction[0] + direction[1] * direction[1];
  const Real lambda_tmp = std::sqrt(norm2 * p.mu * inv_Tr);

  Dv[0] = op(-lambda_tmp);
  Dv[1] = op(lambda_tmp);
  Dv[2] = op(0.0);
}

// ----------------------------------------------------------------------------

template <typename V1, typename M1, typename M2, typename M3>
void AdvDiffLinear2DFOS::flux_jacobian_eigen_structure(const Properties &p, const V1 &direction,
                                                       M1 &Rv, M2 &Lv, M3 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "Rv must be a matrix");
  static_assert(common::IntegersAreEqual<math::TensorRank<M2>::value, math::tensor_rank_2>::value,
                "Lv must be a matrix");
  static_assert(common::IntegersAreEqual<math::TensorRank<M3>::value, math::tensor_rank_2>::value,
                "Dv must be a matrix");

  const Real inv_Tr     = p.mu / (p.Lr * p.Lr);
  const Real norm2      = direction[0] * direction[0] + direction[1] * direction[1];
  const Real norm       = std::sqrt(norm2);
  const Real lambda_tmp = norm * std::sqrt(p.mu * inv_Tr);

  Rv(0, 0) = p.Lr * norm;
  Rv(0, 1) = -p.Lr * norm;
  Rv(0, 2) = 0.0;
  Rv(1, 0) = direction[X0];
  Rv(1, 1) = direction[X0];
  Rv(1, 2) = -direction[X1];
  Rv(2, 0) = direction[X1];
  Rv(2, 1) = direction[X1];
  Rv(2, 2) = direction[X0];

  Lv(0, 0) = 0.5 / (p.Lr * norm);
  Lv(0, 1) = 0.5 * direction[X0] / norm2;
  Lv(0, 2) = 0.5 * direction[X1] / norm2;
  Lv(1, 0) = -0.5 / (p.Lr * norm);
  Lv(1, 1) = 0.5 * direction[X0] / norm2;
  Lv(1, 2) = 0.5 * direction[X1] / norm2;
  Lv(2, 0) = 0.0;
  Lv(2, 1) = -direction[X1] / norm2;
  Lv(2, 2) = direction[X0] / norm2;

  Dv(0, 0) = -lambda_tmp;
  Dv(1, 1) = lambda_tmp;
  Dv(2, 2) = 0.0;
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2, typename M1, typename OP>
void AdvDiffLinear2DFOS::build_K_mat(const Properties &p, const V1 &dir, V2 &eigvalues, M1 &K_mat,
                                     const OP &op)
{
  const Real inv_Tr     = p.mu / (p.Lr * p.Lr);
  const Real norm2      = dir[0] * dir[0] + dir[1] * dir[1];
  const Real norm       = std::sqrt(norm2);
  const Real lambda_tmp = std::sqrt(norm2 * p.mu * inv_Tr);

  const Real L0 = op(-lambda_tmp);
  const Real L1 = op(lambda_tmp);
  const Real L2 = op(0.0);

  eigvalues[0] = L0;
  eigvalues[1] = L1;
  eigvalues[2] = L2;

  K_mat(0, 0) = 0.5 * (L0 + L1);
  K_mat(0, 1) = 0.5 * dir[X0] * p.Lr * (L0 - L1) / norm;
  K_mat(0, 2) = 0.5 * dir[X1] * p.Lr * (L0 - L1) / norm;
  K_mat(1, 0) = 0.5 * dir[X0] * (L0 - L1) / (p.Lr * norm);
  K_mat(1, 1) = 0.5 * (dir[X0] * dir[X0] * (L0 + L1) + 2. * dir[X1] * dir[X1] * L2) / norm2;
  K_mat(1, 2) = 0.5 * (dir[X0] * dir[X1] * (L0 + L1) - 2. * dir[X0] * dir[X1] * L2) / norm2;
  K_mat(2, 0) = 0.5 * dir[X1] * (L0 - L1) / (p.Lr * norm);
  K_mat(2, 1) = 0.5 * (dir[X0] * dir[X1] * (L0 + L1) - 2. * dir[X0] * dir[X1] * L2) / norm2;
  K_mat(2, 2) = 0.5 * (dir[X1] * dir[X1] * (L0 + L1) + 2. * dir[X0] * dir[X0] * L2) / norm2;
}

// ----------------------------------------------------------------------------

template <typename M1, typename V1>
void AdvDiffLinear2DFOS::residual(const Properties &p, std::array<M1, DIM> &flux_jacob, V1 &res)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "M1 must be a matrix type");
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'res' must be a vector");

  M1 &A = flux_jacob[X0];
  M1 &B = flux_jacob[X1];

  const Real Tr = p.Lr * p.Lr / p.mu;

  A(0, 0) = 0.0;
  A(0, 1) = -p.mu;
  A(0, 2) = 0.0;
  A(1, 0) = -1. / Tr;
  A(1, 1) = 0.0;
  A(1, 2) = 0.0;
  A(2, 0) = 0.0;
  A(2, 1) = 0.0;
  A(2, 1) = 0.0;

  B(0, 0) = 0.0;
  B(0, 1) = 0.0;
  B(0, 2) = -p.mu;
  B(1, 0) = 0.0;
  B(1, 1) = 0.0;
  B(1, 2) = 0.0;
  B(2, 0) = -1. / Tr;
  B(2, 1) = 0.0;
  B(2, 1) = 0.0;

  res = A * p.grad_vars.const_col(X) + B * p.grad_vars.const_col(Y);

  // Add source term to the residual
  res[1] += p.vars[1] / Tr;
  res[2] += p.vars[2] / Tr;
}

// ----------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit

#endif
