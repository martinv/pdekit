/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE polynomial_test
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "math/polynomials/Polynomial.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(polynomial_utest_1D)
{
  math::Polynomial<1> linear_poly;
  common::tuple_of<1, Int> lin_exponents;
  std::get<0>(lin_exponents) = 2.0;
  linear_poly.add_term(-2.03, lin_exponents);
  std::cout << linear_poly.to_string() << std::endl;

  common::tuple_of<1, Real> variables;
  std::get<0>(variables) = 4.0;
  std::cout << "Polynomial(" << std::get<0>(variables) << ") = " << linear_poly.evaluate(variables)
            << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(polynomial_utest_2D)
{
  math::Polynomial<2> quad_poly;

  common::tuple_of<2, Int> quad_exponents;

  std::get<0>(quad_exponents) = 2.0;
  std::get<1>(quad_exponents) = 3.0;
  quad_poly.add_term(1.1, quad_exponents);

  std::get<0>(quad_exponents) = 4.0;
  std::get<1>(quad_exponents) = 6.0;
  quad_poly.add_term(2.2, quad_exponents);

  std::get<0>(quad_exponents) = -8.0;
  std::get<1>(quad_exponents) = 12.0;
  quad_poly.add_term(-0.2, quad_exponents);

  std::cout << quad_poly.to_string() << std::endl;

  common::tuple_of<2, Real> variables;
  std::get<0>(variables) = 1.0;
  std::get<1>(variables) = 1.0;
  std::cout << "Polynomial(" << std::get<0>(variables) << "," << std::get<1>(variables)
            << ") = " << quad_poly.evaluate(variables) << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(polynomial_utest_3D)
{
  math::Polynomial<3> cubic_poly;

  common::tuple_of<3, Int> cubic_exponents;

  std::get<0>(cubic_exponents) = 2.0;
  std::get<1>(cubic_exponents) = 3.0;
  std::get<2>(cubic_exponents) = 4.0;
  cubic_poly.add_term(1.1, cubic_exponents);

  std::get<0>(cubic_exponents) = 4.0;
  std::get<1>(cubic_exponents) = 6.0;
  std::get<2>(cubic_exponents) = 8.0;
  cubic_poly.add_term(2.2, cubic_exponents);

  std::get<0>(cubic_exponents) = -4.0;
  std::get<1>(cubic_exponents) = -6.0;
  std::get<1>(cubic_exponents) = -8.0;
  cubic_poly.add_term(-0.2, cubic_exponents);

  std::cout << cubic_poly.to_string() << std::endl;

  common::tuple_of<3, Real> variables;
  std::get<0>(variables) = 1.0;
  std::get<1>(variables) = 1.0;
  std::get<2>(variables) = 2.0;
  std::cout << "Polynomial(" << std::get<0>(variables) << "," << std::get<1>(variables) << ","
            << std::get<2>(variables) << ") = " << cubic_poly.evaluate(variables) << std::endl;

  math::Polynomial<3> cubic_poly2 = 2.0 * cubic_poly;

  std::cout << cubic_poly2.to_string() << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(add_polynomials)
{
  math::Polynomial<3> p1;

  common::tuple_of<3, Int> cubic_exponents;

  std::get<0>(cubic_exponents) = 2.0;
  std::get<1>(cubic_exponents) = 3.0;
  std::get<2>(cubic_exponents) = 4.0;
  p1.add_term(1.1, cubic_exponents);

  std::get<0>(cubic_exponents) = 4.0;
  std::get<1>(cubic_exponents) = 6.0;
  std::get<2>(cubic_exponents) = 8.0;
  p1.add_term(2.2, cubic_exponents);

  std::get<0>(cubic_exponents) = -4.0;
  std::get<1>(cubic_exponents) = -6.0;
  std::get<1>(cubic_exponents) = -8.0;
  p1.add_term(-0.2, cubic_exponents);

  math::Polynomial<3> p2;
  std::get<0>(cubic_exponents) = 2.0;
  std::get<1>(cubic_exponents) = 3.0;
  std::get<2>(cubic_exponents) = 4.0;
  p2.add_term(1.9, cubic_exponents);

  math::Polynomial<3> p3;

  p3 = p1 + p2;

  std::cout << p1.to_string() << std::endl;
  std::cout << p2.to_string() << std::endl;
  std::cout << p3.to_string() << std::endl;
}

// ----------------------------------------------------------------------------
