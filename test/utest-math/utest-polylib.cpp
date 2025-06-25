/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE polylib_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/polynomials/PolyLib.hpp"

using namespace pdekit;

BOOST_AUTO_TEST_CASE(polylib_utest)
{
  const Uint N     = 10;
  const Real alpha = 0.0;
  const Real beta  = 0.0;

  math::DenseDVec<Real> z(N), w(N);

  math::zwgj(z, w, N, alpha, beta);

  std::cout << "z = " << std::endl << z << std::endl;
  std::cout << "w = " << std::endl << w << std::endl;

  Real sum_w = 0.0;
  for (Uint i = 0; i < w.size(); ++i)
  {
    sum_w += w[i];
  }

  std::cout << "Sum of weights = " << sum_w << std::endl;

  std::cout << "Resetting z to zeros ..." << std::endl;
  z.fill(0.0);

  math::jac_zeros(N, z, alpha, beta);

  std::cout << "z recomputed = " << std::endl << z << std::endl;
}
