#include "math/polynomials/JacobiPolynomial.hpp"

#include <iostream>

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

JacobiPolynomial::JacobiPolynomial()
    : a1(0.0), a2(0.0), a3(0.0), a4(0.0), pn_minus1(0.0), pn(0.0), pn_plus1(0.0)
{
}

// ----------------------------------------------------------------------------

JacobiPolynomial::~JacobiPolynomial()
{
}

// ----------------------------------------------------------------------------

Real JacobiPolynomial::operator()(const Uint N, const Real alpha, const Real beta, const Real x)
{

  if (N == 0)
  {
    pn_plus1 = 1.0;
  }
  else if (N == 1)
  {
    pn_plus1 = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
  }
  else
  {
    pn_minus1 = 1.0;
    pn        = 0.5 * (alpha - beta + (alpha + beta + 2.0) * x);
    // The recursive formula gives value for P(n+1), therefore to compute
    // P(2), we have to start with n = 1!
    for (Uint n = 1; n < N; ++n)
    {
      a1 = 2. * (n + 1.) * (n + alpha + beta + 1.) * (2. * n + alpha + beta);
      a2 = (2. * n + alpha + beta + 1.) * (alpha * alpha - beta * beta);
      a3 = (2. * n + alpha + beta) * (2. * n + alpha + beta + 1) * (2. * n + alpha + beta + 2.);
      a4 = 2. * (n + alpha) * (n + beta) * (2. * n + alpha + beta + 2.);

      pn_plus1  = 1.0 / a1 * ((a2 + a3 * x) * pn - a4 * pn_minus1);
      pn_minus1 = pn;
      pn        = pn_plus1;
    }
  }
  return pn_plus1;
}

// ----------------------------------------------------------------------------

Real JacobiPolynomial::dx(const Uint N, const Real alpha, const Real beta, const Real x)
{
  if (N == 0)
  {
    return 0.0;
  }
  else
  {
    operator()(N - 1, alpha + 1, beta + 1, x);
    // return std::sqrt(N*(N+alpha+beta+1))*pn_plus1;

    return 0.5 * (N + alpha + beta + 1.) * pn_plus1;
  }
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit
