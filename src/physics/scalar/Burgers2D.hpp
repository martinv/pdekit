#ifndef PDEKIT_Physics_Burgers2D_hpp
#define PDEKIT_Physics_Burgers2D_hpp

#include <array>

#include "physics/PhysModelT.hpp"
#include "physics/scalar/Adv2DProperties.hpp"

namespace pdekit
{

namespace physics
{

class Burgers2D : public PhysModelT<Adv2DProperties>
{

  public:
  /// TYPEDEFS

  typedef typename PhysModelT<Adv2DProperties>::Properties Properties;
  typedef PhysicsClassScalar physics_class;

  /// Constructor
  Burgers2D();

  /// Destructor
  ~Burgers2D();

  enum
  {
    U = 0
  };

  // ----------------------------------------------------------------------------

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

  // ----------------------------------------------------------------------------

  /// Compute the value of the solution 'u' from the properties
  /// @param p     ... properties
  /// @return vars ... computed variables

  template <typename V1>
  static void compute_variables(const Properties &p, V1 &vars);

  // ----------------------------------------------------------------------------

  /// Compute the fluxes
  /// @param p     ... input physical properties to compute the flux
  /// @return flux ... matrix of fluxes computed from the properties
  template <typename M1>
  static void flux(const Properties &p, M1 &flux);

  // ----------------------------------------------------------------------------

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

  // ----------------------------------------------------------------------------

  /// Compute the eigenvalues of the flux jacobian
  /// @param p         ... physical properties
  /// @param direction ... direction vector, the jacobian that we
  ///                      decompose will be J = direction[X]*dFx/dU +
  /// direction[Y]*dFy/dU
  /// @return Dv       ... vector of computed eigenvalues
  template <typename V1, typename V2>
  static void flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv);

  // ----------------------------------------------------------------------------

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

  // ----------------------------------------------------------------------------

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

  // ----------------------------------------------------------------------------

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

// ----------------------------------------------------------------------------

template <typename V1, typename V2, typename M1>
void Burgers2D::compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars,
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

  p.V[X] = sol[U]; // nonlinear advection velocity in the x-direction
  p.V[Y] = 1.;     //

  p.mu = 0.; // no diffusion
}

// ----------------------------------------------------------------------------

template <typename V1>
void Burgers2D::compute_variables(const Properties &p, V1 &vars)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'vars' must be a vector");
  vars = p.vars;
}

// ----------------------------------------------------------------------------

template <typename M1>
void Burgers2D::flux(const Properties &p, M1 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'flux' must be a matrix");
  flux(0, X) = 0.5 * p.V[X] * p.vars[U];
  flux(0, Y) = p.V[Y] * p.vars[U];
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2>
void Burgers2D::flux(const Properties &p, const V1 &direction, V2 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'flux' must be a vector");

  flux[0] = p.vars[U] * (0.5 * p.V[X] * direction[X] + p.V[Y] * direction[Y]);
}

// --------------------------------------------------------------------------

template <typename V1>
void Burgers2D::source_term(const Properties &p, V1 &src_term)
{
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2>
void Burgers2D::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");
  Dv[0] = p.V[X] * direction[X] + p.V[Y] * direction[Y];
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2, typename OP>
void Burgers2D::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv, OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");
  Dv[0] = op(p.V[X] * direction[X] + p.V[Y] * direction[Y]);
}

// ----------------------------------------------------------------------------

template <typename V1, typename M1, typename M2, typename M3>
void Burgers2D::flux_jacobian_eigen_structure(const Properties &p, const V1 &direction, M1 &Rv,
                                              M2 &Lv, M3 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "Rv must be a matrix");
  static_assert(common::IntegersAreEqual<math::TensorRank<M2>::value, math::tensor_rank_2>::value,
                "Lv must be a matrix");
  static_assert(common::IntegersAreEqual<math::TensorRank<M3>::value, math::tensor_rank_2>::value,
                "Dv must be a matrix");

  Rv(0, 0) = 1.;
  Lv(0, 0) = 1.;
  Dv(0, 0) = p.V[X] * direction[X] + p.V[Y] * direction[Y];
}

// ----------------------------------------------------------------------------

template <typename V1, typename V2, typename M1, typename OP>
void Burgers2D::build_K_mat(const Properties &p, const V1 &direction, V2 &eigvalues, M1 &K_mat,
                            const OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'eigvalues' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'K_mat' must be a matrix");

  eigvalues[0] = p.V[X] * direction[X] + p.V[Y] * direction[Y];
  K_mat(0, 0)  = 1. * op(p.V[X] * direction[X] + p.V[Y] * direction[Y]) * 1.;
}

// --------------------------------------------------------------------------

template <typename M1, typename V1>
void Burgers2D::residual(const Properties &p, std::array<M1, DIM> &flux_jacob, V1 &res)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "M1 must be a matrix type");
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'res' must be a vector");

  M1 &A = flux_jacob[X];
  M1 &B = flux_jacob[Y];

  A(0, 0) = p.V[X];
  B(0, 0) = p.V[Y];

  //   lambda_1 * du/dx        + lambda_2 * du/dy
  //      y     * du/dx        +   (-x)   * du/dy
  res = A * p.grad_vars.const_col(X) + B * p.grad_vars.const_col(Y);
}

// ----------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit

#endif
