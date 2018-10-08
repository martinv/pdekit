#ifndef PDEKIT_Physics_Euler_3D_Cons_hpp
#define PDEKIT_Physics_Euler_3D_Cons_hpp

#include <array>

#include "physics/PhysModelT.hpp"
#include "physics/euler/Euler3DProperties.hpp"

namespace pdekit
{

namespace physics
{

// ----------------------------------------------------------------------------

class Euler3DCons : public PhysModelT<Euler3DProperties>
{
  public:
  /// TYPEDEFS

  typedef typename PhysModelT<Euler3DProperties>::Properties Properties;
  typedef PhysicsClassEuler physics_class;

  /// Constructor
  Euler3DCons();

  /// Destructor
  ~Euler3DCons();

  enum
  {
    Rho  = 0,
    RhoU = 1,
    RhoV = 2,
    RhoW = 3,
    RhoE = 4
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
void Euler3DCons::compute_properties(const V1 &coord, const V2 &sol, const M1 &grad_vars,
                                     Properties &p)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'coord' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'sol' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'grad_vars' must be a matrix");

  p.coords    = coord;     // cache the coordiantes locally
  p.vars      = sol;       // cache the variables locally
  p.grad_vars = grad_vars; // cache the gradient of variables locally

  p.rho  = sol[Rho];
  p.rhou = sol[RhoU];
  p.rhov = sol[RhoV];
  p.rhow = sol[RhoW];
  p.rhoE = sol[RhoE];

  p.inv_rho = 1. / p.rho;

  p.V[X] = p.rhou * p.inv_rho;
  p.V[Y] = p.rhov * p.inv_rho;
  p.V[Z] = p.rhow * p.inv_rho;

  p.v2 = p.V[X] * p.V[X] + p.V[Y] * p.V[Y] + p.V[Z] * p.V[Z];

  p.P = p.gamma_minus_1 * (p.rhoE - 0.5 * p.rho * p.v2);

  if (p.P <= 0.)
  {
    std::cout << "rho     : " << p.rho << std::endl;
    std::cout << "rhou    : " << p.rhou << std::endl;
    std::cout << "rhov    : " << p.rhov << std::endl;
    std::cout << "rhow    : " << p.rhow << std::endl;
    std::cout << "rhoE    : " << p.rhoE << std::endl;
    std::cout << "P       : " << p.P << std::endl;
    std::cout << "u       : " << p.V[X] << std::endl;
    std::cout << "v       : " << p.V[Y] << std::endl;
    std::cout << "w       : " << p.V[Z] << std::endl;
    std::cout << "v2  : " << p.v2 << std::endl;

    std::cout << "Pressure is negative at coordinates [" << coord[X] << ",";
    std::cout << coord[Y] << "," << coord[Z] << "]" << std::endl;
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
void Euler3DCons::compute_variables(const Properties &p, V1 &vars)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'vars' must be a vector");
  vars[Rho]  = p.rho;
  vars[RhoU] = p.rhou;
  vars[RhoV] = p.rhov;
  vars[RhoW] = p.rhow;
  vars[RhoE] = p.rhoE;
}

// --------------------------------------------------------------------------

template <typename M1>
void Euler3DCons::flux(const Properties &p, M1 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'flux' must be a matrix");
  flux(0, X) = p.rhou;                // rho.u
  flux(1, X) = p.rhou * p.V[X] + p.P; // rho.u^2 + P
  flux(2, X) = p.rhou * p.V[Y];       // rho.u.v
  flux(3, X) = p.rhou * p.V[Z];       // rho.u.w
  flux(4, X) = p.rhou * p.H;          // rho.u.H

  flux(0, Y) = p.rhov;                // rho.v
  flux(1, Y) = p.rhov * p.V[X];       // rho.v.u
  flux(2, Y) = p.rhov * p.V[Y] + p.P; // rho.v^2 + P
  flux(3, Y) = p.rhov * p.V[Z];       // rho.v.w
  flux(4, Y) = p.rhov * p.H;          // rho.v.H

  flux(0, Z) = p.rhow;                // rho.w
  flux(1, Z) = p.rhow * p.V[X];       // rho.w.u
  flux(2, Z) = p.rhow * p.V[Y];       // rho.w.v
  flux(3, Z) = p.rhow * p.V[Z] + p.P; // rho.w^2 + P
  flux(4, Z) = p.rhow * p.H;          // rho.w.H
}

// --------------------------------------------------------------------------

template <typename V1, typename V2>
void Euler3DCons::flux(const Properties &p, const V1 &direction, V2 &flux)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'flux' must be a vector");

  const Real rhoum = p.rhou * direction[X] + p.rhov * direction[Y] + p.rhow * direction[Z];

  flux[0] = rhoum;
  flux[1] = rhoum * p.V[X] + p.P * direction[X];
  flux[2] = rhoum * p.V[Y] + p.P * direction[Y];
  flux[3] = rhoum * p.V[Z] + p.P * direction[Z];
  flux[4] = rhoum * p.H;
}

// --------------------------------------------------------------------------

template <typename V1>
inline void Euler3DCons::source_term(const Properties &p, V1 &src_term)
{
}

// --------------------------------------------------------------------------

template <typename V1, typename V2>
void Euler3DCons::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");
  const Real um = p.V[X] * direction[X] + p.V[Y] * direction[Y] + p.V[Z] * direction[Z];

  Dv[0] = um;
  Dv[1] = um;
  Dv[2] = um;
  Dv[3] = um + p.a;
  Dv[4] = um - p.a;
}

// --------------------------------------------------------------------------

template <typename V1, typename V2, typename OP>
void Euler3DCons::flux_jacobian_eigen_values(const Properties &p, const V1 &direction, V2 &Dv,
                                             OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'Dv' must be a vector");
  const Real um = p.V[X] * direction[X] + p.V[Y] * direction[Y] + p.V[Z] * direction[Z];

  const Real op_um = op(um);

  Dv[0] = op_um;
  Dv[1] = op_um;
  Dv[2] = op_um;
  Dv[3] = op_um + p.a;
  Dv[4] = op_um - p.a;
}

// --------------------------------------------------------------------------

template <typename V1, typename M1, typename M2, typename M3>
void Euler3DCons::flux_jacobian_eigen_structure(const Properties &p, const V1 &direction, M1 &Rv,
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

  const Real norm = std::sqrt(direction[X] * direction[X] + direction[Y] * direction[Y] +
                              direction[Z] * direction[Z]);
  //  if (norm < 1.e-9)
  //  {
  //    std::cerr << "Error, zero norm: " << norm << std::endl;
  //    std::cin.get();
  //  }

  const Real nx = (norm == 0) ? direction[X] : direction[X] / norm;
  const Real ny = (norm == 0) ? direction[Y] : direction[Y] / norm;
  const Real nz = (norm == 0) ? direction[Z] : direction[Z] / norm;

  const Real inv_a  = 1. / p.a;
  const Real inv_a2 = inv_a * inv_a;

  const Real um = p.V[X] * nx + p.V[Y] * ny + p.V[Z] * nz;
  const Real ra = 0.5 * p.rho * inv_a;

  const Real gu_a  = p.gamma_minus_1 * p.V[X] * inv_a;
  const Real gv_a  = p.gamma_minus_1 * p.V[Y] * inv_a;
  const Real gw_a  = p.gamma_minus_1 * p.V[Z] * inv_a;
  const Real rho_a = p.rho * p.a;

  const Real gm1_ov_rhoa = p.gamma_minus_1 / rho_a;

  const Real ke = 0.5 * p.v2;
  // const Real gm2 = p.half_gm1_v2 * inv_a2;
  const Real gm2 = (0.5 * p.gamma_minus_1 * p.v2) * inv_a2;
  const Real k2  = 1.0 - gm2;
  const Real k3  = -p.gamma_minus_1 * inv_a2;

  // matrix of right eigen vectors R

  Rv(0, 0) = nx;
  Rv(0, 1) = ny;
  Rv(0, 2) = nz;
  Rv(0, 3) = ra;
  Rv(0, 4) = ra;

  Rv(1, 0) = p.V[X] * nx;
  Rv(1, 1) = p.V[X] * ny - p.rho * nz;
  Rv(1, 2) = p.V[X] * nz + p.rho * ny;
  Rv(1, 3) = ra * (p.V[X] + p.a * nx);
  Rv(1, 4) = ra * (p.V[X] - p.a * nx);

  Rv(2, 0) = p.V[Y] * nx + p.rho * nz;
  Rv(2, 1) = p.V[Y] * ny;
  Rv(2, 2) = p.V[Y] * nz - p.rho * nx;
  Rv(2, 3) = ra * (p.V[Y] + p.a * ny);
  Rv(2, 4) = ra * (p.V[Y] - p.a * ny);

  Rv(3, 0) = p.V[Z] * nx - p.rho * ny;
  Rv(3, 1) = p.V[Z] * ny + p.rho * nx;
  Rv(3, 2) = p.V[Z] * nz;
  Rv(3, 3) = ra * (p.V[Z] + p.a * nz);
  Rv(3, 4) = ra * (p.V[Z] - p.a * nz);

  Rv(4, 0) = ke * nx + p.rho * (p.V[Y] * nz - p.V[Z] * ny);
  Rv(4, 1) = ke * ny + p.rho * (p.V[Z] * nx - p.V[X] * nz);
  Rv(4, 2) = ke * nz + p.rho * (p.V[X] * ny - p.V[Y] * nx);
  Rv(4, 3) = ra * (p.H + p.a * um);
  Rv(4, 4) = ra * (p.H - p.a * um);

  // matrix of left eigen vectors L = R.inverse();

  Lv(0, 0) = nx * k2 - p.inv_rho * (p.V[Y] * nz - p.V[Z] * ny);
  Lv(0, 1) = gu_a * inv_a * nx;
  Lv(0, 2) = gv_a * inv_a * nx + nz * p.inv_rho;
  Lv(0, 3) = gw_a * inv_a * nx - ny * p.inv_rho;
  Lv(0, 4) = k3 * nx;

  Lv(1, 0) = ny * k2 - p.inv_rho * (p.V[Z] * nx - p.V[X] * nz);
  Lv(1, 1) = gu_a * inv_a * ny - nz * p.inv_rho;
  Lv(1, 2) = gv_a * inv_a * ny;
  Lv(1, 3) = gw_a * inv_a * ny + nx * p.inv_rho;
  Lv(1, 4) = k3 * ny;

  Lv(2, 0) = nz * k2 - p.inv_rho * (p.V[X] * ny - p.V[Y] * nx);
  Lv(2, 1) = gu_a * inv_a * nz + ny * p.inv_rho;
  Lv(2, 2) = gv_a * inv_a * nz - nx * p.inv_rho;
  Lv(2, 3) = gw_a * inv_a * nz;
  Lv(2, 4) = k3 * nz;

  Lv(3, 0) = p.a * p.inv_rho * (gm2 - um / p.a);
  Lv(3, 1) = p.inv_rho * (nx - gu_a);
  Lv(3, 2) = p.inv_rho * (ny - gv_a);
  Lv(3, 3) = p.inv_rho * (nz - gw_a);
  Lv(3, 4) = gm1_ov_rhoa;

  Lv(4, 0) = p.a * p.inv_rho * (gm2 + um / p.a);
  Lv(4, 1) = p.inv_rho * (-nx - gu_a);
  Lv(4, 2) = p.inv_rho * (-ny - gv_a);
  Lv(4, 3) = p.inv_rho * (-nz - gw_a);
  Lv(4, 4) = gm1_ov_rhoa;

  // diagonal matrix of eigenvalues

  Dv(0, 0) = norm * um;
  Dv(1, 1) = norm * um;
  Dv(2, 2) = norm * um;
  Dv(3, 3) = norm * (um + p.a);
  Dv(4, 4) = norm * (um - p.a);
}

// --------------------------------------------------------------------------

template <typename V1, typename V2, typename M1, typename OP>
void Euler3DCons::build_K_mat(const Properties &p, const V1 &direction, V2 &eigvalues, M1 &K_mat,
                              const OP &op)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'direction' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<V2>::value, math::tensor_rank_1>::value,
                "'eigvalues' must be a vector");
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "'K_mat' must be a matrix");

  const Real norm = std::sqrt(direction[X] * direction[X] + direction[Y] * direction[Y] +
                              direction[Z] * direction[Z]);
  //  if (norm < 1.e-9)
  //  {
  //    std::cerr << "Error, zero norm: " << norm << std::endl;
  //    std::cin.get();
  //  }

  const Real nx = (norm == 0) ? direction[X] : direction[X] / norm;
  const Real ny = (norm == 0) ? direction[Y] : direction[Y] / norm;
  const Real nz = (norm == 0) ? direction[Z] : direction[Z] / norm;

  const Real inv_a  = 1. / p.a;
  const Real inv_a2 = inv_a * inv_a;

  const Real um = p.V[X] * nx + p.V[Y] * ny + p.V[Z] * nz;
  const Real ra = 0.5 * p.rho * inv_a;

  const Real gu_a  = p.gamma_minus_1 * p.V[X] * inv_a;
  const Real gv_a  = p.gamma_minus_1 * p.V[Y] * inv_a;
  const Real gw_a  = p.gamma_minus_1 * p.V[Z] * inv_a;
  const Real rho_a = p.rho * p.a;

  const Real gm1_ov_rhoa = p.gamma_minus_1 / rho_a;

  const Real ke = 0.5 * p.v2;
  // const Real gm2 = p.half_gm1_v2 * inv_a2;
  const Real gm2 = (0.5 * p.gamma_minus_1 * p.v2) * inv_a2;
  const Real k2  = 1.0 - gm2;
  const Real k3  = -p.gamma_minus_1 * inv_a2;

  eigvalues[0] = norm * um;
  eigvalues[1] = norm * um;
  eigvalues[2] = norm * um;
  eigvalues[3] = norm * (um + p.a);
  eigvalues[4] = norm * (um - p.a);

  const Real lambda0 = op(eigvalues[0]);
  const Real lambda1 = op(eigvalues[1]);
  const Real lambda2 = op(eigvalues[2]);
  const Real lambda3 = op(eigvalues[3]);
  const Real lambda4 = op(eigvalues[4]);

  // matrix of right eigen vectors R
  const Real Rv00 = nx;
  const Real Rv01 = ny;
  const Real Rv02 = nz;
  const Real Rv03 = ra;
  const Real Rv04 = ra;

  const Real Rv10 = p.V[X] * nx;
  const Real Rv11 = p.V[X] * ny - p.rho * nz;
  const Real Rv12 = p.V[X] * nz + p.rho * ny;
  const Real Rv13 = ra * (p.V[X] + p.a * nx);
  const Real Rv14 = ra * (p.V[X] - p.a * nx);

  const Real Rv20 = p.V[Y] * nx + p.rho * nz;
  const Real Rv21 = p.V[Y] * ny;
  const Real Rv22 = p.V[Y] * nz - p.rho * nx;
  const Real Rv23 = ra * (p.V[Y] + p.a * ny);
  const Real Rv24 = ra * (p.V[Y] - p.a * ny);

  const Real Rv30 = p.V[Z] * nx - p.rho * ny;
  const Real Rv31 = p.V[Z] * ny + p.rho * nx;
  const Real Rv32 = p.V[Z] * nz;
  const Real Rv33 = ra * (p.V[Z] + p.a * nz);
  const Real Rv34 = ra * (p.V[Z] - p.a * nz);

  const Real Rv40 = ke * nx + p.rho * (p.V[Y] * nz - p.V[Z] * ny);
  const Real Rv41 = ke * ny + p.rho * (p.V[Z] * nx - p.V[X] * nz);
  const Real Rv42 = ke * nz + p.rho * (p.V[X] * ny - p.V[Y] * nx);
  const Real Rv43 = ra * (p.H + p.a * um);
  const Real Rv44 = ra * (p.H - p.a * um);

  // matrix of left eigen vectors L = R.inverse();

  const Real Lv00 = nx * k2 - p.inv_rho * (p.V[Y] * nz - p.V[Z] * ny);
  const Real Lv01 = gu_a * inv_a * nx;
  const Real Lv02 = gv_a * inv_a * nx + nz * p.inv_rho;
  const Real Lv03 = gw_a * inv_a * nx - ny * p.inv_rho;
  const Real Lv04 = k3 * nx;

  const Real Lv10 = ny * k2 - p.inv_rho * (p.V[Z] * nx - p.V[X] * nz);
  const Real Lv11 = gu_a * inv_a * ny - nz * p.inv_rho;
  const Real Lv12 = gv_a * inv_a * ny;
  const Real Lv13 = gw_a * inv_a * ny + nx * p.inv_rho;
  const Real Lv14 = k3 * ny;

  const Real Lv20 = nz * k2 - p.inv_rho * (p.V[X] * ny - p.V[Y] * nx);
  const Real Lv21 = gu_a * inv_a * nz + ny * p.inv_rho;
  const Real Lv22 = gv_a * inv_a * nz - nx * p.inv_rho;
  const Real Lv23 = gw_a * inv_a * nz;
  const Real Lv24 = k3 * nz;

  const Real Lv30 = p.a * p.inv_rho * (gm2 - um / p.a);
  const Real Lv31 = p.inv_rho * (nx - gu_a);
  const Real Lv32 = p.inv_rho * (ny - gv_a);
  const Real Lv33 = p.inv_rho * (nz - gw_a);
  const Real Lv34 = gm1_ov_rhoa;

  const Real Lv40 = p.a * p.inv_rho * (gm2 + um / p.a);
  const Real Lv41 = p.inv_rho * (-nx - gu_a);
  const Real Lv42 = p.inv_rho * (-ny - gv_a);
  const Real Lv43 = p.inv_rho * (-nz - gw_a);
  const Real Lv44 = gm1_ov_rhoa;

  K_mat(0, 0) = Rv00 * lambda0 * Lv00 + Rv01 * lambda1 * Lv10 + Rv02 * lambda2 * Lv20 +
                Rv03 * lambda3 * Lv30 + Rv04 * lambda4 * Lv40;

  K_mat(0, 1) = Rv00 * lambda0 * Lv01 + Rv01 * lambda1 * Lv11 + Rv02 * lambda2 * Lv21 +
                Rv03 * lambda3 * Lv31 + Rv04 * lambda4 * Lv41;

  K_mat(0, 2) = Rv00 * lambda0 * Lv02 + Rv01 * lambda1 * Lv12 + Rv02 * lambda2 * Lv22 +
                Rv03 * lambda3 * Lv32 + Rv04 * lambda4 * Lv42;

  K_mat(0, 3) = Rv00 * lambda0 * Lv03 + Rv01 * lambda1 * Lv13 + Rv02 * lambda2 * Lv23 +
                Rv03 * lambda3 * Lv33 + Rv04 * lambda4 * Lv43;

  K_mat(0, 4) = Rv00 * lambda0 * Lv04 + Rv01 * lambda1 * Lv14 + Rv02 * lambda2 * Lv24 +
                Rv03 * lambda3 * Lv34 + Rv04 * lambda4 * Lv44;

  // -----------------------------

  K_mat(1, 0) = Rv10 * lambda0 * Lv00 + Rv11 * lambda1 * Lv10 + Rv12 * lambda2 * Lv20 +
                Rv13 * lambda3 * Lv30 + Rv14 * lambda4 * Lv40;

  K_mat(1, 1) = Rv10 * lambda0 * Lv01 + Rv11 * lambda1 * Lv11 + Rv12 * lambda2 * Lv21 +
                Rv13 * lambda3 * Lv31 + Rv14 * lambda4 * Lv41;

  K_mat(1, 2) = Rv10 * lambda0 * Lv02 + Rv11 * lambda1 * Lv12 + Rv12 * lambda2 * Lv22 +
                Rv13 * lambda3 * Lv32 + Rv14 * lambda4 * Lv42;

  K_mat(1, 3) = Rv10 * lambda0 * Lv03 + Rv11 * lambda1 * Lv13 + Rv12 * lambda2 * Lv23 +
                Rv13 * lambda3 * Lv33 + Rv14 * lambda4 * Lv43;

  K_mat(1, 4) = Rv10 * lambda0 * Lv04 + Rv11 * lambda1 * Lv14 + Rv12 * lambda2 * Lv24 +
                Rv13 * lambda3 * Lv34 + Rv14 * lambda4 * Lv44;

  // -----------------------------

  K_mat(2, 0) = Rv20 * lambda0 * Lv00 + Rv21 * lambda1 * Lv10 + Rv22 * lambda2 * Lv20 +
                Rv23 * lambda3 * Lv30 + Rv24 * lambda4 * Lv40;

  K_mat(2, 1) = Rv20 * lambda0 * Lv01 + Rv21 * lambda1 * Lv11 + Rv22 * lambda2 * Lv21 +
                Rv23 * lambda3 * Lv31 + Rv24 * lambda4 * Lv41;

  K_mat(2, 2) = Rv20 * lambda0 * Lv02 + Rv21 * lambda1 * Lv12 + Rv22 * lambda2 * Lv22 +
                Rv23 * lambda3 * Lv32 + Rv24 * lambda4 * Lv42;

  K_mat(2, 3) = Rv20 * lambda0 * Lv03 + Rv21 * lambda1 * Lv13 + Rv22 * lambda2 * Lv23 +
                Rv23 * lambda3 * Lv33 + Rv24 * lambda4 * Lv43;

  K_mat(2, 4) = Rv20 * lambda0 * Lv04 + Rv21 * lambda1 * Lv14 + Rv22 * lambda2 * Lv24 +
                Rv23 * lambda3 * Lv34 + Rv24 * lambda4 * Lv44;

  // -----------------------------

  K_mat(3, 0) = Rv30 * lambda0 * Lv00 + Rv31 * lambda1 * Lv10 + Rv32 * lambda2 * Lv20 +
                Rv33 * lambda3 * Lv30 + Rv34 * lambda4 * Lv40;

  K_mat(3, 1) = Rv30 * lambda0 * Lv01 + Rv31 * lambda1 * Lv11 + Rv32 * lambda2 * Lv21 +
                Rv33 * lambda3 * Lv31 + Rv34 * lambda4 * Lv41;

  K_mat(3, 2) = Rv30 * lambda0 * Lv02 + Rv31 * lambda1 * Lv12 + Rv32 * lambda2 * Lv22 +
                Rv33 * lambda3 * Lv32 + Rv34 * lambda4 * Lv42;

  K_mat(3, 3) = Rv30 * lambda0 * Lv03 + Rv31 * lambda1 * Lv13 + Rv32 * lambda2 * Lv23 +
                Rv33 * lambda3 * Lv33 + Rv34 * lambda4 * Lv43;

  K_mat(3, 4) = Rv30 * lambda0 * Lv04 + Rv31 * lambda1 * Lv14 + Rv32 * lambda2 * Lv24 +
                Rv33 * lambda3 * Lv34 + Rv34 * lambda4 * Lv44;

  // -----------------------------

  K_mat(4, 0) = Rv40 * lambda0 * Lv00 + Rv41 * lambda1 * Lv10 + Rv42 * lambda2 * Lv20 +
                Rv43 * lambda3 * Lv30 + Rv44 * lambda4 * Lv40;

  K_mat(4, 1) = Rv40 * lambda0 * Lv01 + Rv41 * lambda1 * Lv11 + Rv42 * lambda2 * Lv21 +
                Rv43 * lambda3 * Lv31 + Rv44 * lambda4 * Lv41;

  K_mat(4, 2) = Rv40 * lambda0 * Lv02 + Rv41 * lambda1 * Lv12 + Rv42 * lambda2 * Lv22 +
                Rv43 * lambda3 * Lv32 + Rv44 * lambda4 * Lv42;

  K_mat(4, 3) = Rv40 * lambda0 * Lv03 + Rv41 * lambda1 * Lv13 + Rv42 * lambda2 * Lv23 +
                Rv43 * lambda3 * Lv33 + Rv44 * lambda4 * Lv43;

  K_mat(4, 4) = Rv40 * lambda0 * Lv04 + Rv41 * lambda1 * Lv14 + Rv42 * lambda2 * Lv24 +
                Rv43 * lambda3 * Lv34 + Rv44 * lambda4 * Lv44;
}

// --------------------------------------------------------------------------

template <typename M1, typename V1>
void Euler3DCons::residual(const Properties &p, std::array<M1, DIM> &flux_jacob, V1 &res)
{
  static_assert(common::IntegersAreEqual<math::TensorRank<M1>::value, math::tensor_rank_2>::value,
                "M1 must be a matrix type");
  static_assert(common::IntegersAreEqual<math::TensorRank<V1>::value, math::tensor_rank_1>::value,
                "'res' must be a vector");

  const Real gamma_minus_3 = p.gamma - 3.;
  const Real half_gm1_v2   = 0.5 * p.gamma_minus_1 * p.v2;

  const Real uu = p.V[X] * p.V[X];
  const Real uv = p.V[X] * p.V[Y];
  const Real uw = p.V[X] * p.V[Z];
  const Real vv = p.V[Y] * p.V[Y];
  const Real vw = p.V[Y] * p.V[Z];
  const Real ww = p.V[Z] * p.V[Z];

  JM &Jx = flux_jacob[X];

  Jx.fill(0.0);

  //    A.setZero(); // assume are zeroed

  //  Jx(0,0) = 0.0;
  Jx(0, 1) = 1.0;
  //  Jx(0,2) = 0.0;
  //  Jx(0,3) = 0.0;
  //  Jx(0,4) = 0.0;
  Jx(1, 0) = half_gm1_v2 - uu;
  Jx(1, 1) = -gamma_minus_3 * p.V[X];
  Jx(1, 2) = -p.gamma_minus_1 * p.V[Y];
  Jx(1, 3) = -p.gamma_minus_1 * p.V[Z];
  Jx(1, 4) = p.gamma_minus_1;
  Jx(2, 0) = -uv;
  Jx(2, 1) = p.V[Y];
  Jx(2, 2) = p.V[X];
  //  Jx(2,3) = 0.0;
  //  Jx(2,4) = 0.0;
  Jx(3, 0) = -uw;
  Jx(3, 1) = p.V[Z];
  //  Jx(3,2) = 0.0;
  Jx(3, 3) = p.V[X];
  //  Jx(3,4) = 0.0;
  Jx(4, 0) = p.V[X] * (half_gm1_v2 - p.H);
  Jx(4, 1) = -p.gamma_minus_1 * uu + p.H;
  Jx(4, 2) = -p.gamma_minus_1 * uv;
  Jx(4, 3) = -p.gamma_minus_1 * uw;
  Jx(4, 4) = p.gamma * p.V[X];

  JM &Jy = flux_jacob[Y];

  Jy.fill(0.0);

  //    B.setZero(); // assume are zeroed

  //  Jy(0,0) = 0.0;
  //  Jy(0,1) = 0.0;
  Jy(0, 2) = 1.0;
  //  Jy(0,3) = 0.0;
  //  Jy(0,4) = 0.0;
  Jy(1, 0) = -uv;
  Jy(1, 1) = p.V[Y];
  Jy(1, 2) = p.V[X];
  //  Jy(1,3) = 0.0;
  //  Jy(1,4) = 0.0;
  Jy(2, 0) = half_gm1_v2 - vv;
  Jy(2, 1) = -p.gamma_minus_1 * p.V[X];
  Jy(2, 2) = -gamma_minus_3 * p.V[Y];
  Jy(2, 3) = -p.gamma_minus_1 * p.V[Z];
  Jy(2, 4) = p.gamma_minus_1;
  Jy(3, 0) = -vw;
  //  Jy(3,1) = 0.0;
  Jy(3, 2) = p.V[Z];
  Jy(3, 3) = p.V[Y];
  //  Jy(3,4) = 0.0;
  Jy(4, 0) = p.V[Y] * (half_gm1_v2 - p.H);
  Jy(4, 1) = -p.gamma_minus_1 * uv;
  Jy(4, 2) = -p.gamma_minus_1 * vv + p.H;
  Jy(4, 3) = -p.gamma_minus_1 * vw;
  Jy(4, 4) = p.gamma * p.V[Y];

  JM &Jz = flux_jacob[Z];

  Jz.fill(0.0);

  //  Jz(0,0) = 0.0;
  //  Jz(0,1) = 0.0;
  //  Jz(0,2) = 0.0;
  Jz(0, 3) = 1.0;
  //  Jz(0,4) = 0.0;
  Jz(1, 0) = -uw;
  Jz(1, 1) = p.V[Z];
  //  Jz(1,2) = 0.0;
  Jz(1, 3) = p.V[X];
  //  Jz(1,4) = 0.0;
  Jz(2, 0) = -vw;
  //  Jz(2,1) = 0.0;
  Jz(2, 2) = p.V[Z];
  Jz(2, 3) = p.V[Y];
  //  Jz(2,4) = 0.0;
  Jz(3, 0) = half_gm1_v2 - ww;
  Jz(3, 1) = -p.gamma_minus_1 * p.V[X];
  Jz(3, 2) = -p.gamma_minus_1 * p.V[Y];
  Jz(3, 3) = -gamma_minus_3 * p.V[Z];
  Jz(3, 4) = p.gamma_minus_1;
  Jz(4, 0) = p.V[Z] * (half_gm1_v2 - p.H);
  Jz(4, 1) = -p.gamma_minus_1 * uw;
  Jz(4, 2) = -p.gamma_minus_1 * vw;
  Jz(4, 3) = -p.gamma_minus_1 * ww + p.H;
  Jz(4, 4) = p.gamma * p.V[Z];

  res =
      Jx * p.grad_vars.const_col(X) + Jy * p.grad_vars.const_col(Y) + Jz * p.grad_vars.const_col(Z);
}

// --------------------------------------------------------------------------

} // namespace physics

} // namespace pdekit

#endif
