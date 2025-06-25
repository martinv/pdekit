#include "math/polynomials/Polynomial.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

namespace poly_detail
{

// ----------------------------------------------------------------------------
// Specialized printing class for 1D polynomials
// ----------------------------------------------------------------------------

std::string PolyToString<1>::apply(const std::vector<Real> &coeffs,
                                   const std::vector<Int> &exponents)
{
  std::stringstream ss;

  const Uint n_terms = coeffs.size();

  if ((n_terms > 0) && (coeffs[0] != 0.0))
  {
    ss << coeffs[0] << "*x^" << exponents[0];
  }

  for (Uint t = 1; t < n_terms; ++t)
  {
    if (coeffs[t] != 0.0)
    {
      if (coeffs[t] > 0.0)
      {
        ss << " + " << coeffs[t];
        ;
      }
      else
      {
        ss << " - " << std::abs(coeffs[t]);
      }
      ss << coeffs[t] << "*x^" << exponents[t];
    }
  }

  std::string result(ss.str());
  return result;
}

// ----------------------------------------------------------------------------
// Specialized printing class for 2D polynomials
// ----------------------------------------------------------------------------

std::string PolyToString<2>::apply(const std::vector<Real> &coeffs,
                                   const std::vector<Int> &exponents)
{
  std::stringstream ss;
  const char var_names[3] = {'x', 'y', 'z'};

  const Uint n_terms = coeffs.size();

  if ((n_terms > 0) && (coeffs[0] != 0.0))
  {
    ss << coeffs[0] << "*x^" << exponents[0] << "*y^" << exponents[1];
  }

  for (Uint t = 1; t < n_terms; ++t)
  {
    if (coeffs[t] != 0.0)
    {
      if (coeffs[t] > 0.0)
      {
        ss << " + " << coeffs[t];
      }
      else
      {
        ss << " - " << std::abs(coeffs[t]);
      }
      for (Uint v = 0; v < 2; ++v)
      {
        ss << "*" << var_names[v] << "^";
        if (exponents[2 * t + v] < 0.0)
        {
          ss << "(" << exponents[2 * t + v] << ")";
        }
        else
        {
          ss << exponents[2 * t + v];
        }
      }
    }
  }

  std::string result(ss.str());
  return result;
}

// ----------------------------------------------------------------------------
// Specialized printing class for 3D polynomials
// ----------------------------------------------------------------------------

std::string PolyToString<3>::apply(const std::vector<Real> &coeffs,
                                   const std::vector<Int> &exponents)
{
  std::stringstream ss;
  const char var_names[3] = {'x', 'y', 'z'};

  const Uint n_terms = coeffs.size();

  if ((n_terms > 0) && (coeffs[0] != 0.0))
  {
    ss << coeffs[0] << "*x^" << exponents[0] << "*y^" << exponents[1] << "*z^" << exponents[2];
  }

  for (Uint t = 1; t < n_terms; ++t)
  {
    if (coeffs[t] != 0.0)
    {
      if (coeffs[t] > 0.0)
      {
        ss << " + " << coeffs[t];
      }
      else
      {
        ss << " - " << std::abs(coeffs[t]);
      }

      for (Uint v = 0; v < 3; ++v)
      {
        ss << "*" << var_names[v] << "^";
        if (exponents[3 * t + v] < 0.0)
        {
          ss << "(" << exponents[3 * t + v] << ")";
        }
        else
        {
          ss << exponents[3 * t + v];
        }
      }
    }
  }

  std::string result(ss.str());
  return result;
}

// ----------------------------------------------------------------------------

} // namespace poly_detail

} // namespace math

} // namespace pdekit
