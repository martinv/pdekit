#ifndef PDEKIT_Physics_Scalar_Navier_Stokes_2D_Properties_hpp
#define PDEKIT_Physics_Scalar_Navier_Stokes_2D_Properties_hpp

#include "common/Constants.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace physics
{

struct Euler2DProperties
{
  /// TYPEDEFS:

  enum
  {
    DIM = _2D
  }; /// number of dimensions
  enum
  {
    NEQ = 4
  }; /// number of independent variables or equations

  typedef math::DenseSVec<Real, DIM> CoordV;        // type of geometry coordinates
                                                    // vector
  typedef math::DenseDMat<Real> SolM;               // type of solution matrix (should be
                                                    // sized to nb_nodes x NEQ)
  typedef math::DenseSVec<Real, NEQ> SolV;          // type of solution variables vector
  typedef math::DenseSVec<Real, NEQ> FluxV;         // type of flux vector
  typedef math::DenseSVec<Real, NEQ> SolGradV;      // type of solution gradient
                                                    // vector
  typedef math::DenseSMat<Real, NEQ, DIM> SolGradM; // type of solution
                                                    // gradient matrix

  /// Constructor
  Euler2DProperties();

  /// Destructor
  ~Euler2DProperties();

  /// Print the contents:
  void print() const;

  CoordV coords;      /// position in domain
  SolV vars;          /// independent variables with positions described in coords
  SolGradM grad_vars; /// gradient of independent variables

  Real gamma;         /// specific heat ratio
  Real R;             /// gas constant
  Real gamma_minus_1; /// specific heat ratio minus one, very commonly used

  Real rho;  /// density
  Real rhou; /// rho.u
  Real rhov; /// rho.v
  Real rhoE; /// rho.E

  Real inv_rho; /// inverse of density, very commonly used

  CoordV V; /// velocity
  Real v2;  /// u^2 + v^2

  Real H; /// specific enthalpy
  // Real a2;          /// square of speed of sound, very commonly used
  Real a; /// speed of sound
  Real P; /// pressure
  Real T; /// temperature
  Real E; /// specific internal energy
  // Real half_gm1_v2; /// 1/2.(g-1).(u^2+v^2), very commonly used
  Real Ma; /// mach number
};

} // namespace physics

} // namespace pdekit

#endif
