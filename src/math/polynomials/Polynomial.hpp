#ifndef PDEKIT_Math_Polynomial_hpp
#define PDEKIT_Math_Polynomial_hpp

#include "common/PDEKit.hpp"
#include "common/TupleOf.hpp"
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

namespace poly_detail
{
template <size_t NVars>
class PolyToString
{
  public:
  static std::string apply(const std::vector<Real> &coeffs, const std::vector<Int> &exponents)
  {
    std::stringstream ss;

    const Uint n_terms = coeffs.size();

    if ((n_terms > 0) && (coeffs[0] != 0.0))
    {
      ss << coeffs[0];
      for (Uint v = 0; v < NVars; ++v)
      {
        ss << "*x" << v << "^" << exponents[v];
      }
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
        // ss << m_coeffs[t];
        for (Uint v = 0; v < NVars; ++v)
        {
          ss << "*x" << v << "^" << exponents[NVars * t + v];
        }
      }
    }

    std::string result(ss.str());
    return result;
  }
};

// ----------------------------------------------------------------------------

template <>
class PolyToString<1>
{
  public:
  static std::string apply(const std::vector<Real> &coeffs, const std::vector<Int> &exponents);
};

template <>
class PolyToString<2>
{
  public:
  static std::string apply(const std::vector<Real> &coeffs, const std::vector<Int> &exponents);
};

template <>
class PolyToString<3>
{
  public:
  static std::string apply(const std::vector<Real> &coeffs, const std::vector<Int> &exponents);
};

// ----------------------------------------------------------------------------

} // namespace poly_detail

// ----------------------------------------------------------------------------
// A class representing polynomials with real coefficients
// and integer exponents
// ----------------------------------------------------------------------------

template <size_t NVars>
class Polynomial;

template <size_t NVars>
Polynomial<NVars> operator*(const Real coeff, const Polynomial<NVars> &poly);

template <size_t NVars>
Polynomial<NVars> operator+(const Polynomial<NVars> &lhs, const Polynomial<NVars> &rhs);

// ----------------------------------------------------------------------------

template <size_t NVars>
class Polynomial
{
  public:
  /// Default constructor
  Polynomial() = default;

  /// Copy constructor
  Polynomial(const Polynomial &rhs);

  /// Move constructor
  Polynomial(const Polynomial &&rhs);

  /// Default destructor
  ~Polynomial() = default;

  /// Assignment operator
  Polynomial &operator=(const Polynomial &rhs);

  /// Move assignment operator
  Polynomial &operator=(Polynomial &&rhs);

  /// Add a new term
  void add_term(const Real coeff, const common::tuple_of<NVars, Real> &exponents);

  /// Evaluate polynomial for given values of variables
  Real evaluate(const common::tuple_of<NVars, Real> &variables) const;

  friend Polynomial operator*<NVars>(const Real coeff, const Polynomial &poly);

  friend Polynomial operator+<NVars>(const Polynomial &lhs, const Polynomial &rhs);

  /// Generate a string representation of the polynomial
  std::string to_string() const;

  private:
  /// Helper classes to process tuples of exponents on input
  template <size_t NArgs, size_t ArgPos>
  struct ProcessPolyRepresentation
  {
    static void read_exponents(const common::tuple_of<NArgs, Int> &in_exponents,
                               std::vector<Int> &out_exponents)
    {
      out_exponents.push_back(std::get<ArgPos>(in_exponents));
      ProcessPolyRepresentation<NArgs, ArgPos + 1>::read_exponents(in_exponents, out_exponents);
    }

    static void evaluate_term(const common::tuple_of<NArgs, Real> &vars,
                              const std::vector<Int> &exponents, const Uint term_pos, Real &term)
    {
      term *= std::pow(std::get<ArgPos>(vars), exponents[NArgs * term_pos + ArgPos]);
      ProcessPolyRepresentation<NArgs, ArgPos + 1>::evaluate_term(vars, exponents, term_pos, term);
    }
  };

  template <size_t NArgs>
  struct ProcessPolyRepresentation<NArgs, NArgs - 1>
  {
    static void read_exponents(const common::tuple_of<NArgs, Int> &in_exponents,
                               std::vector<Int> &out_exponents)
    {
      out_exponents.push_back(std::get<NArgs - 1>(in_exponents));
    }

    static void evaluate_term(const common::tuple_of<NArgs, Real> &vars,
                              const std::vector<Int> &exponents, const Uint term_pos, Real &term)
    {
      term *= std::pow(std::get<NArgs - 1>(vars), exponents[NArgs * term_pos + NArgs - 1]);
    }
  };

  /// String conversion classes
  template <size_t NumberVars, size_t Dummy>
  class PolyToString;

  /// DATA
  /// Coefficient of each term in polynomial
  std::vector<Real> m_coeffs;

  /// Exponents for all variables of each term, grouped by terms
  /// Example: if this is a polynomial with 2 variables (x and y),
  /// then the exponents are stored as pairs { [n0x, n0y], [n1x, n1y],
  /// [n2x, n2y], ..., [nkx, nky] }
  /// The polynomial is then
  /// m_coeffs[0] * x^n0x * y^n0y + m_coeffs[1] * x^n1x * y^n1y + ...
  /// + m_coeffs[k] * x^nkx * y^nky
  std::vector<Int> m_exponents;
};

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars>::Polynomial(const Polynomial &rhs)
{
  m_coeffs.resize(rhs.m_coeffs.size());
  m_exponents.resize(rhs.m_exponents.size());

  m_coeffs    = rhs.m_coeffs;
  m_exponents = rhs.m_exponents;
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars>::Polynomial(const Polynomial &&rhs)
    : m_coeffs(std::move(rhs.m_coeffs)), m_exponents(std::move(rhs.m_exponents))
{
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars> &Polynomial<NVars>::operator=(const Polynomial &rhs)
{
  m_coeffs.resize(rhs.m_coeffs.size());
  m_exponents.resize(rhs.m_exponents.size());

  m_coeffs    = rhs.m_coeffs;
  m_exponents = rhs.m_exponents;

  return *this;
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars> &Polynomial<NVars>::operator=(Polynomial &&rhs)
{
  m_coeffs    = std::move(rhs.m_coeffs);
  m_exponents = std::move(rhs.m_exponents);

  return *this;
}

// ----------------------------------------------------------------------------

template <size_t NVars>
void Polynomial<NVars>::add_term(const Real coeff, const common::tuple_of<NVars, Real> &exponents)
{
  m_coeffs.push_back(coeff);
  ProcessPolyRepresentation<NVars, 0>::read_exponents(exponents, m_exponents);
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Real Polynomial<NVars>::evaluate(const common::tuple_of<NVars, Real> &variables) const
{
  Real result = 0;
  for (Uint t = 0; t < m_coeffs.size(); ++t)
  {
    Real term = m_coeffs[t];
    ProcessPolyRepresentation<NVars, 0>::evaluate_term(variables, m_exponents, t, term);
    result += term;
  }
  return result;
}

// ----------------------------------------------------------------------------

template <size_t NVars>
std::string Polynomial<NVars>::to_string() const
{
  return poly_detail::PolyToString<NVars>::apply(m_coeffs, m_exponents);
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars> operator*(const Real coeff, const Polynomial<NVars> &poly)
{
  Polynomial<NVars> result(poly);

  for (Uint i = 0; i < result.m_coeffs.size(); ++i)
  {
    result.m_coeffs[i] *= coeff;
  }

  return result;
}

// ----------------------------------------------------------------------------

template <size_t NVars>
Polynomial<NVars> operator+(const Polynomial<NVars> &lhs, const Polynomial<NVars> &rhs)
{
  Polynomial<NVars> result(lhs);

  // Now add add coefficients and exponents from RHS
  for (Uint i = 0; i < rhs.m_coeffs.size(); ++i)
  {
    Uint matching_term_pos = result.m_coeffs.size() + 1;
    for (Uint j = 0; j < result.m_coeffs.size(); ++j)
    {
      bool exponents_equal = true;

      for (Uint k = 0; k < NVars; ++k)
      {
        if (result.m_exponents[j * NVars + k] != rhs.m_exponents[i * NVars + k])
        {
          exponents_equal = false;
          break;
        }
      }

      if (exponents_equal)
      {
        matching_term_pos = j;
        break;
      }

    } // loop over terms of polynomial 'result'

    if (matching_term_pos < result.m_coeffs.size())
    {
      result.m_coeffs[matching_term_pos] += rhs.m_coeffs[i];
    }
    else
    {
      result.m_coeffs.push_back(rhs.m_coeffs[i]);
      for (Uint k = 0; k < NVars; ++k)
      {
        result.m_exponents.push_back(rhs.m_exponents[i * NVars + k]);
      }
    } // else
  }   // Loop over terms of rhs polynomial

  return result;
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
