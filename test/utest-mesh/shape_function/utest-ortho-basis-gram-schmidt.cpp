/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ortho_basis_gram_schmidt_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "math/algo/GramSchmidt.hpp"
#include "math/polynomials/Polynomial.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

struct OrthoBasisGramSchmidt_Fixture
{
  /// common setup for each test case
  OrthoBasisGramSchmidt_Fixture();

  /// common tear-down for each test case
  ~OrthoBasisGramSchmidt_Fixture();
};

OrthoBasisGramSchmidt_Fixture::OrthoBasisGramSchmidt_Fixture()
{
}

OrthoBasisGramSchmidt_Fixture::~OrthoBasisGramSchmidt_Fixture()
{
}

// ----------------------------------------------------------------------------

class ScalarProductPoly1D
{
  public:
  Real operator()(const math::Polynomial<1> &p1, const math::Polynomial<1> &p2) const
  {
    const Real xq[3] = {-0.7745966692414833770358531, 0.0, 0.7745966692414833770358531};
    const Real wq[3] = {5. / 9., 8. / 9., 5. / 9};

    common::tuple_of<1, Real> variables;
    Real result = 0.0;

    for (Uint q = 0; q < 3; ++q)
    {
      std::get<0>(variables) = xq[q];
      result += wq[q] * p1.evaluate(variables) * p2.evaluate(variables);
    }
    return result;
  }
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(OrthoBasisGramSchmidt_TestSuite, OrthoBasisGramSchmidt_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(basis_gram_schmidt_1D_utest)
{
  // ScalarProductPoly1D scal_prod;
  std::vector<math::Polynomial<1>> poly_in;
  std::vector<math::Polynomial<1>> poly_out;

  poly_in.resize(3);
  poly_out.resize(3);

  common::tuple_of<1, Int> p2_exponents;

  std::get<0>(p2_exponents) = 0.0;
  poly_in[0].add_term(1.0, p2_exponents);

  std::get<0>(p2_exponents) = 1.0;
  poly_in[1].add_term(1.0, p2_exponents);

  std::get<0>(p2_exponents) = 2.0;
  poly_in[2].add_term(1.0, p2_exponents);

  std::cout << "Basis 0 = " << poly_in[0].to_string() << std::endl;
  std::cout << "Basis 1 = " << poly_in[1].to_string() << std::endl;
  std::cout << "Basis 2 = " << poly_in[2].to_string() << std::endl;

  // math::algo::GramSchmidt(scal_prod, poly_in, poly_out);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
