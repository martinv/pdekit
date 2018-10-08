#ifndef PDEKIT_Physics_Scalar_Adv2D_Properties_hpp
#define PDEKIT_Physics_Scalar_Adv2D_Properties_hpp

#include "common/Constants.hpp"
#include "math/DenseConstVecView.hpp"
#include "math/DenseDMat.hpp"
#include "math/DenseSMat.hpp"
#include "math/DenseSVec.hpp"

namespace pdekit
{

namespace physics
{

struct Adv2DProperties
{
  /// TYPEDEFS:

  enum
  {
    DIM = _2D
  }; /// number of dimensions
  enum
  {
    NEQ = 1
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
  Adv2DProperties();

  /// Destructor
  ~Adv2DProperties();

  /// Print the contents:
  void print() const;

  CoordV coords;      /// position in domain
  SolV vars;          /// independent variables with positions described in coords
  SolGradM grad_vars; /// gradient of independent variables
  CoordV V;           /// advection speed
  Real mu;            /// scalar diffusion coefficient
};

} // namespace physics

} // namespace pdekit

#endif
