#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "common/PDEKit.hpp"
#include "common/StringUtils.hpp"
#include "math/polynomials/PolyLib.hpp"

using namespace pdekit;

class LagrangeInterpolant
{
  public:
  LagrangeInterpolant(const Uint order, const std::string point_distribution);

  ~LagrangeInterpolant();

  Real value(const Real x) const;

  Uint order() const;

  private:
  const Uint m_order;
  math::DenseDVec<Real> m_pts;
};

// ============================================================================

LagrangeInterpolant::LagrangeInterpolant(const Uint order, const std::string point_distribution)
    : m_order(order)
{
  m_pts.resize(order + 1);

  if (point_distribution == "Equidist")
  {
    const Real dx = 2. / order;
    for (Uint i = 0; i < (order + 1); ++i)
    {
      m_pts[i] = -1. + i * dx;
    }
  }
  else if (point_distribution == "Chebyshev")
  {
    // alpha = beta = -0.5 => Chebyshev polynomial as special case
    // of Jacobi polynomial
    math::jac_zeros(m_order + 1, m_pts, -0.5, -0.5);
  }
}

// ============================================================================

LagrangeInterpolant::~LagrangeInterpolant()
{
}

// ============================================================================

Real LagrangeInterpolant::value(const Real x) const
{
  Real value = 0.0;

  for (Uint i = 0; i < m_pts.size(); ++i)
  {
    Real denominator = 1.0;
    Real nominator   = 1.0;

    for (Uint j = 0; j < m_pts.size(); ++j)
    {
      if (i != j)
      {
        denominator *= (m_pts[i] - m_pts[j]);
        nominator *= (x - m_pts[j]);
      }
    }

    const Real Lagr_i = nominator / denominator;

    value += 1. / (1. + 25. * m_pts[i] * m_pts[i]) * Lagr_i;
  }
  return value;
}

// ============================================================================

Uint LagrangeInterpolant::order() const
{
  return m_order;
}

// ============================================================================

int main()
{
  const Uint P = 10;

  LagrangeInterpolant L1(P, "Equidist");
  LagrangeInterpolant L2(P, "Chebyshev");

  const Uint Npts = 101;
  const Real dx   = 2.0 / (Npts - 1);

  std::ofstream outfile;
  outfile.precision(10);
  outfile.setf(std::ios::fixed);

  outfile.open("Runge_phenomenon_P" + common::StringUtils::to_string(P) + ".dat");

  for (Uint i = 0; i < Npts; ++i)
  {
    const Real x = -1. + i * dx;

    const Real exact = 1. / (1. + 25 * x * x);

    outfile << std::setw(14) << x << " " << std::setw(14) << exact << " " << std::setw(14)
            << L1.value(x) << " " << std::setw(14) << L2.value(x) << std::endl;
  }

  outfile.close();

  return 0;
}
