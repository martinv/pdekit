#ifndef PDEKIT_Physics_Euler_2D_Cons_hpp
#define PDEKIT_Physics_Euler_2D_Cons_hpp

#include <array>

#include "physics/PhysModelT.hpp"
#include "physics/euler/Euler2DProperties.hpp"

namespace pdekit
{

namespace physics
{

// --------------------------------------------------------------------------====

class Euler2DCons : public PhysModelT<Euler2DProperties>
{
  public:
  /// TYPEDEFS

  typedef typename PhysModelT<Euler2DProperties>::Properties Properties;
  typedef PhysicsClassEuler physics_class;

  /// Constructor
  Euler2DCons();

  /// Destructor
  ~Euler2DCons();

  enum
  {
    Rho  = 0,
    RhoU = 1,
    RhoV = 2,
    RhoE = 3
  };

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
  ///                      direction[Y]*dFy/dU
  /// @return Dv       ... vector of computed eigenvalues
  /// @param  Op       ... extra operator applied to the eigen values (max,
  /// '+' operator
  ///                      for example)
  template <typename V1, typename V2, typename OP>
  static void flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv, OP &op);

  // --------------------------------------------------------------------------

  /// Decompose the eigen structure of the flux jacobians projected on the
  /// gradients
  /// @param p         ... physical properties
  /// @param direction ... direction onto which the Jacobian
  ///                      is projected (J = direction[X]*dFx/dU +
  ///                      direction[Y]*dFy/dU)
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
void Euler2DCons::compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars,
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

  p.rho  = sol[Rho];
  p.rhou = sol[RhoU];
  p.rhov = sol[RhoV];
  p.rhoE = sol[RhoE];

  p.inv_rho = 1. / p.rho;

  p.V[X] = p.rhou * p.inv_rho;
  p.V[Y] = p.rhov * p.inv_rho;

  p.v2 = p.V[X] * p.V[X] + p.V[Y] * p.V[Y];

  p.P = p.gamma_minus_1 * (p.rhoE - 0.5 * p.rho * p.v2);

  if (p.P <= 0.)
  {
    std::cout << "rho   : " << p.rho << std::endl;
    std::cout << "rhou  : " << p.rhou << std::endl;
    std::cout << "rhov  : " << p.rhov << std::endl;
    std::cout << "rhoE  : " << p.rhoE << std::endl;
    std::cout << "P     : " << p.P << std::endl;
    std::cout << "u     : " << p.V[X] << std::endl;
    std::cout << "v     : " << p.V[Y] << std::endl;
    std::cout << "uuvv  : " << p.v2 << std::endl;

    std::cout << "Pressure is negative at coordinates [" << coord[X] << ",";
    std::cout << coord[Y] << "]" << std::endl;
  }

  const Real RT = p.P * p.inv_rho; // RT = p/rho

  p.E = p.rhoE * p.inv_rho; // E = rhoE / rho

  p.H = p.E + RT; // H = E + p/rho

  p.a = sqrt(p.gamma * RT);

  // p.a2 = p.a * p.a;

  p.Ma = sqrt(p.v2) / p.a;

  p.T = RT / p.R;

  // p.half_gm1_v2 = 0.5 * p.gamma_minus_1 * p.v2;
}

// --------------------------------------------------------------------------

template <typename V1>
void Euler2DCons::compute_variables(const Properties &p, V1 &vars)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'vars' must be a vector");
  vars[Rho]  = p.rho;
  vars[RhoU] = p.rhou;
  vars[RhoV] = p.rhov;
  vars[RhoE] = p.rhoE;
}

// --------------------------------------------------------------------------

template <typename M1>
void Euler2DCons::flux(const Properties &p, M1 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'flux' must be a matrix");
  flux(0, X) = p.rhou;                // rho.u
  flux(1, X) = p.rhou * p.V[X] + p.P; // rho.u^2 + P
  flux(2, X) = p.rhou * p.V[Y];       // rho.u.v
  flux(3, X) = p.rhou * p.H;          // rho.u.H

  flux(0, Y) = p.rhov;                // rho.v
  flux(1, Y) = p.rhov * p.V[X];       // rho.v.u
  flux(2, Y) = p.rhov * p.V[Y] + p.P; // rho.v^2 + P
  flux(3, Y) = p.rhov * p.H;          // rho.v.H
}

// --------------------------------------------------------------------------

template <typename V1, typename V2>
void Euler2DCons::flux(const Properties &p, const V1 &direction, V2 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'flux' must be a vector");

  const Real rhoum = p.rhou * direction[X] + p.rhov * direction[Y];

  flux[0] = rhoum;
  flux[1] = rhoum * p.V[X] + p.P * direction[X];
  flux[2] = rhoum * p.V[Y] + p.P * direction[Y];
  flux[3] = rhoum * p.H;
}

// --------------------------------------------------------------------------

template <typename V1>
inline void Euler2DCons::source_term(const Properties &p, V1 &src_term)
{
}

// --------------------------------------------------------------------------

template <typename V1, typename V2>
void Euler2DCons::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");

  const Real um = p.V[X] * direction[X] + p.V[Y] * direction[Y];

  Dv[0] = um;
  Dv[1] = um;
  Dv[2] = um + p.a;
  Dv[3] = um - p.a;
}

// --------------------------------------------------------------------------

template <typename V1, typename V2, typename OP>
void Euler2DCons::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv,
                                             OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");
  const Real um = p.V[X] * direction[X] + p.V[Y] * direction[Y];

  const Real op_um = op(um);

  Dv[0] = op_um;
  Dv[1] = op_um;
  Dv[2] = op_um + p.a;
  Dv[3] = op_um - p.a;
}

// --------------------------------------------------------------------------

template <typename V1, typename M1, typename M2, typename M3>
void Euler2DCons::flux_jacobian_eigen_structure(const Properties &p, const V1 &direction, M1 &Rv,
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

  const Real norm = std::sqrt(direction[X] * direction[X] + direction[Y] * direction[Y]);
  //  if (norm < 1.e-9)
  //  {
  //    std::cerr << "Error, zero norm: " << norm << std::endl;
  //    std::cin.get();
  //  }

  const Real nx = (norm < 1.e-9) ? direction[X] : direction[X] / norm;
  const Real ny = (norm < 1.e-9) ? direction[Y] : direction[Y] / norm;

  const Real inv_a  = 1. / p.a;
  const Real inv_a2 = inv_a * inv_a;

  const Real um = p.V[X] * nx + p.V[Y] * ny;
  const Real ra = 0.5 * p.rho * inv_a;

  // const Real coeffM2 = p.half_gm1_v2 * inv_a2;
  const Real coeffM2 = (0.5 * (p.gamma - 1.0) * p.v2) * inv_a2;

  const Real uDivA = p.gamma_minus_1 * p.V[X] * inv_a;
  const Real vDivA = p.gamma_minus_1 * p.V[Y] * inv_a;
  const Real rho_a = p.rho * p.a;

  const Real gm1_ov_rhoa = p.gamma_minus_1 / rho_a;

  // matrix of right eigenvectors R

  Rv(0, 0) = 1.;
  Rv(0, 1) = 0.;
  Rv(0, 2) = ra;
  Rv(0, 3) = ra;
  Rv(1, 0) = p.V[X];
  Rv(1, 1) = p.rho * ny;
  Rv(1, 2) = ra * (p.V[X] + p.a * nx);
  Rv(1, 3) = ra * (p.V[X] - p.a * nx);
  Rv(2, 0) = p.V[Y];
  Rv(2, 1) = -p.rho * nx;
  Rv(2, 2) = ra * (p.V[Y] + p.a * ny);
  Rv(2, 3) = ra * (p.V[Y] - p.a * ny);
  Rv(3, 0) = 0.5 * p.v2;
  Rv(3, 1) = p.rho * (p.V[X] * ny - p.V[Y] * nx);
  Rv(3, 2) = ra * (p.H + p.a * um);
  Rv(3, 3) = ra * (p.H - p.a * um);

  // matrix of left eigenvectors L = R.inverse();

  Lv(0, 0) = 1. - coeffM2;
  Lv(0, 1) = uDivA * inv_a;
  Lv(0, 2) = vDivA * inv_a;
  Lv(0, 3) = -p.gamma_minus_1 * inv_a2;
  Lv(1, 0) = p.inv_rho * (p.V[Y] * nx - p.V[X] * ny);
  Lv(1, 1) = p.inv_rho * ny;
  Lv(1, 2) = -p.inv_rho * nx;
  Lv(1, 3) = 0.0;
  Lv(2, 0) = p.a * p.inv_rho * (coeffM2 - um * inv_a);
  Lv(2, 1) = p.inv_rho * (nx - uDivA);
  Lv(2, 2) = p.inv_rho * (ny - vDivA);
  Lv(2, 3) = gm1_ov_rhoa;
  Lv(3, 0) = p.a * p.inv_rho * (coeffM2 + um * inv_a);
  Lv(3, 1) = -p.inv_rho * (nx + uDivA);
  Lv(3, 2) = -p.inv_rho * (ny + vDivA);
  Lv(3, 3) = gm1_ov_rhoa;

  // diagonal matrix of eigenvalues

  Dv(0, 0) = norm * um;
  Dv(1, 1) = norm * um;
  Dv(2, 2) = norm * (um + p.a);
  Dv(3, 3) = norm * (um - p.a);
}

// --------------------------------------------------------------------------

template <typename V1, typename V2, typename M1, typename OP>
void Euler2DCons::build_K_mat(const Properties &p, const V1 &direction, V2 &eigvalues, M1 &K_mat,
                              const OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'eigvalues' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'K_mat' must be a matrix");

  const Real norm = std::sqrt(direction[X] * direction[X] + direction[Y] * direction[Y]);
  //  if (norm < 1.e-9)
  //  {
  //    std::cerr << "Error, zero norm: " << norm << std::endl;
  //    std::cin.get();
  //  }

  const Real nx = (norm < 1.e-9) ? direction[X] : direction[X] / norm;
  const Real ny = (norm < 1.e-9) ? direction[Y] : direction[Y] / norm;

  const Real inv_a  = 1. / p.a;
  const Real inv_a2 = inv_a * inv_a;

  const Real um = p.V[X] * nx + p.V[Y] * ny;
  const Real ra = 0.5 * p.rho * inv_a;

  // const Real coeffM2 = p.half_gm1_v2 * inv_a2;
  const Real coeffM2 = (0.5 * (p.gamma - 1.0) * p.v2) * inv_a2;
  const Real uDivA   = p.gamma_minus_1 * p.V[X] * inv_a;
  const Real vDivA   = p.gamma_minus_1 * p.V[Y] * inv_a;
  const Real rho_a   = p.rho * p.a;

  const Real gm1_ov_rhoa = p.gamma_minus_1 / rho_a;

  eigvalues[0] = norm * um;
  eigvalues[1] = norm * um;
  eigvalues[2] = norm * (um + p.a);
  eigvalues[3] = norm * (um - p.a);

  const Real lambda0 = op(eigvalues[0]);
  const Real lambda1 = op(eigvalues[1]);
  const Real lambda2 = op(eigvalues[2]);
  const Real lambda3 = op(eigvalues[3]);

  // matrix of right eigenvectors R
  const Real Rv00 = 1.;
  const Real Rv01 = 0.;
  const Real Rv02 = ra;
  const Real Rv03 = ra;
  const Real Rv10 = p.V[X];
  const Real Rv11 = p.rho * ny;
  const Real Rv12 = ra * (p.V[X] + p.a * nx);
  const Real Rv13 = ra * (p.V[X] - p.a * nx);
  const Real Rv20 = p.V[Y];
  const Real Rv21 = -p.rho * nx;
  const Real Rv22 = ra * (p.V[Y] + p.a * ny);
  const Real Rv23 = ra * (p.V[Y] - p.a * ny);
  const Real Rv30 = 0.5 * p.v2;
  const Real Rv31 = p.rho * (p.V[X] * ny - p.V[Y] * nx);
  const Real Rv32 = ra * (p.H + p.a * um);
  const Real Rv33 = ra * (p.H - p.a * um);

  // matrix of left eigenvectors L = R.inverse();
  const Real Lv00 = 1. - coeffM2;
  const Real Lv01 = uDivA * inv_a;
  const Real Lv02 = vDivA * inv_a;
  const Real Lv03 = -p.gamma_minus_1 * inv_a2;
  const Real Lv10 = p.inv_rho * (p.V[Y] * nx - p.V[X] * ny);
  const Real Lv11 = p.inv_rho * ny;
  const Real Lv12 = -p.inv_rho * nx;
  const Real Lv13 = 0.0;
  const Real Lv20 = p.a * p.inv_rho * (coeffM2 - um * inv_a);
  const Real Lv21 = p.inv_rho * (nx - uDivA);
  const Real Lv22 = p.inv_rho * (ny - vDivA);
  const Real Lv23 = gm1_ov_rhoa;
  const Real Lv30 = p.a * p.inv_rho * (coeffM2 + um * inv_a);
  const Real Lv31 = -p.inv_rho * (nx + uDivA);
  const Real Lv32 = -p.inv_rho * (ny + vDivA);
  const Real Lv33 = gm1_ov_rhoa;

  K_mat(0, 0) =
      Rv00 * lambda0 * Lv00 + Rv01 * lambda1 * Lv10 + Rv02 * lambda2 * Lv20 + Rv03 * lambda3 * Lv30;

  K_mat(0, 1) =
      Rv00 * lambda0 * Lv01 + Rv01 * lambda1 * Lv11 + Rv02 * lambda2 * Lv21 + Rv03 * lambda3 * Lv31;

  K_mat(0, 2) =
      Rv00 * lambda0 * Lv02 + Rv01 * lambda1 * Lv12 + Rv02 * lambda2 * Lv22 + Rv03 * lambda3 * Lv32;

  K_mat(0, 3) =
      Rv00 * lambda0 * Lv03 + Rv01 * lambda1 * Lv13 + Rv02 * lambda2 * Lv23 + Rv03 * lambda3 * Lv33;

  // -----------------------------

  K_mat(1, 0) =
      Rv10 * lambda0 * Lv00 + Rv11 * lambda1 * Lv10 + Rv12 * lambda2 * Lv20 + Rv13 * lambda3 * Lv30;

  K_mat(1, 1) =
      Rv10 * lambda0 * Lv01 + Rv11 * lambda1 * Lv11 + Rv12 * lambda2 * Lv21 + Rv13 * lambda3 * Lv31;

  K_mat(1, 2) =
      Rv10 * lambda0 * Lv02 + Rv11 * lambda1 * Lv12 + Rv12 * lambda2 * Lv22 + Rv13 * lambda3 * Lv32;

  K_mat(1, 3) =
      Rv10 * lambda0 * Lv03 + Rv11 * lambda1 * Lv13 + Rv12 * lambda2 * Lv23 + Rv13 * lambda3 * Lv33;

  // -----------------------------

  K_mat(2, 0) =
      Rv20 * lambda0 * Lv00 + Rv21 * lambda1 * Lv10 + Rv22 * lambda2 * Lv20 + Rv23 * lambda3 * Lv30;

  K_mat(2, 1) =
      Rv20 * lambda0 * Lv01 + Rv21 * lambda1 * Lv11 + Rv22 * lambda2 * Lv21 + Rv23 * lambda3 * Lv31;

  K_mat(2, 2) =
      Rv20 * lambda0 * Lv02 + Rv21 * lambda1 * Lv12 + Rv22 * lambda2 * Lv22 + Rv23 * lambda3 * Lv32;

  K_mat(2, 3) =
      Rv20 * lambda0 * Lv03 + Rv21 * lambda1 * Lv13 + Rv22 * lambda2 * Lv23 + Rv23 * lambda3 * Lv33;

  // -----------------------------

  K_mat(3, 0) =
      Rv30 * lambda0 * Lv00 + Rv31 * lambda1 * Lv10 + Rv32 * lambda2 * Lv20 + Rv33 * lambda3 * Lv30;

  K_mat(3, 1) =
      Rv30 * lambda0 * Lv01 + Rv31 * lambda1 * Lv11 + Rv32 * lambda2 * Lv21 + Rv33 * lambda3 * Lv31;

  K_mat(3, 2) =
      Rv30 * lambda0 * Lv02 + Rv31 * lambda1 * Lv12 + Rv32 * lambda2 * Lv22 + Rv33 * lambda3 * Lv32;

  K_mat(3, 3) =
      Rv30 * lambda0 * Lv03 + Rv31 * lambda1 * Lv13 + Rv32 * lambda2 * Lv23 + Rv33 * lambda3 * Lv33;
}

// --------------------------------------------------------------------------

template <typename M1, typename V1>
void Euler2DCons::residual(const Properties &p, std::array<M1, DIM> &flux_jacob, V1 &res)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "M1 must be a matrix type");
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'res' must be a vector");

  const Real gamma_minus_3 = p.gamma - 3.;
  const Real half_gm1_v2   = 0.5 * (p.gamma - 1.0) * p.v2;

  const Real uu = p.V[X] * p.V[X];
  const Real uv = p.V[X] * p.V[Y];
  const Real vv = p.V[Y] * p.V[Y];

  JM &Jx = flux_jacob[X];

  //    A.setZero(); // assume are zeroed

  Jx(0, 0) = 0.;
  Jx(0, 1) = 1.;
  Jx(0, 2) = 0.;
  Jx(0, 3) = 0.;

  Jx(1, 0) = half_gm1_v2 - uu;
  Jx(1, 1) = -gamma_minus_3 * p.V[X];
  Jx(1, 2) = -p.gamma_minus_1 * p.V[Y];
  Jx(1, 3) = p.gamma_minus_1;

  Jx(2, 0) = -uv;
  Jx(2, 1) = p.V[Y];
  Jx(2, 2) = p.V[X];
  Jx(2, 3) = 0.;

  Jx(3, 0) = half_gm1_v2 * p.V[X] - p.V[X] * p.H;
  Jx(3, 1) = -p.gamma_minus_1 * uu + p.H;
  Jx(3, 2) = -p.gamma_minus_1 * uv;
  Jx(3, 3) = p.gamma * p.V[X];

  JM &Jy = flux_jacob[Y];

  //    B.setZero(); // assume are zeroed

  Jy(0, 0) = 0.;
  Jy(0, 1) = 0.;
  Jy(0, 2) = 1.;
  Jy(0, 3) = 0.;

  Jy(1, 0) = -uv;
  Jy(1, 1) = p.V[Y];
  Jy(1, 2) = p.V[X];
  Jy(1, 3) = 0.;

  Jy(2, 0) = half_gm1_v2 - vv;
  Jy(2, 1) = -p.gamma_minus_1 * p.V[X];
  Jy(2, 2) = -gamma_minus_3 * p.V[Y];
  Jy(2, 3) = p.gamma_minus_1;

  Jy(3, 0) = half_gm1_v2 * p.V[Y] - p.V[Y] * p.H;
  Jy(3, 1) = -p.gamma_minus_1 * uv;
  Jy(3, 2) = -p.gamma_minus_1 * vv + p.H;
  Jy(3, 3) = p.gamma * p.V[Y];

  res = Jx * p.grad_vars.const_col(X) + Jy * p.grad_vars.const_col(Y);
}

// --------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit

#endif
