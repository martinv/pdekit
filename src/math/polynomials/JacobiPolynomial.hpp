#ifndef PDEKIT_Math_Jacobi_Polynomial_hpp
#define PDEKIT_Math_Jacobi_Polynomial_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

/// Class representing Jacobi polynomial (see the book of Sherwin and
/// Karniadakis for definition)
/// Note that for alpha = beta = 0, we have Legendre polynomial
/// For alpha = beta = -1/2, we have Chebyshev polynomial

class JacobiPolynomial
{

  public:
  /// Default constructor
  JacobiPolynomial();

  /// Destructor
  ~JacobiPolynomial();

  /// Evaluate Jacobi polynomial J(alpha,beta,N) at point x
  Real operator()(const Uint N, const Real alpha, const Real beta, const Real x);

  /// Evaluate the derivative of J(alpha,beta,N) at x
  Real dx(const Uint N, const Real alpha, const Real beta, const Real x);

  private:
  /// Helper variables
  Real a1, a2, a3, a4;
  Real pn_minus1;
  Real pn;
  Real pn_plus1;
};

} // namespace math

} // namespace pdekit

#endif
