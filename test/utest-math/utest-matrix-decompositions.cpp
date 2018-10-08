/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE matrix_decompositions_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/DenseDMat.hpp"
#include "math/DenseSMat.hpp"
#include "math/decompositions/EigenvalueDecomposition.hpp"
#include "math/decompositions/LUDecomposition.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

struct MatrixDecompositionsUtestFixture
{
  enum
  {
    rank = 3
  };

  // Fixture constructor
  MatrixDecompositionsUtestFixture();

  // Fixture destructor
  ~MatrixDecompositionsUtestFixture();

  math::DenseDMat<Real> dA, dEig, dR, dRinv, dCheck, dL, dU;
  math::DenseSMat<Real, rank, rank> sA, sEig, sR, sRinv, sCheck, sL, sU;
};

MatrixDecompositionsUtestFixture::MatrixDecompositionsUtestFixture()
{
  dA.resize(rank, rank);
  dA = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(3.);

  dEig.resize(rank, rank);
  dR.resize(rank, rank);
  dRinv.resize(rank, rank);
  dL.resize(rank, rank);
  dU.resize(rank, rank);
  dCheck.resize(rank, rank);

  dEig.fill(0.0);
  dR.fill(0.0);
  dRinv.fill(0.0);
  dL.fill(0.0);
  dU.fill(0.0);
  dCheck.fill(0.0);

  sA = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(3.);
  sEig.fill(0.0);
  sR.fill(0.0);
  sRinv.fill(0.0);
  sL.fill(0.0);
  sU.fill(0.0);
  sCheck.fill(0.0);
}

MatrixDecompositionsUtestFixture::~MatrixDecompositionsUtestFixture()
{
}

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MatrixDecompositions_TestSuite, MatrixDecompositionsUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(eigenvalue_decomposition_utest)
{
  math::EigenvalueDecomposition<Real> eig_decomp;

  // eig_decomp.set_rank(rank);

  eig_decomp.factorize(dA, dEig, dR, dRinv);
  dCheck = dR * dEig * dRinv;

  const Real tol = 1.e-12;
  BOOST_CHECK_CLOSE(dCheck(0, 0), 1.0, tol);
  BOOST_CHECK_CLOSE(dCheck(0, 1), 2.0, tol);
  BOOST_CHECK_CLOSE(dCheck(0, 2), 3.0, tol);
  BOOST_CHECK_CLOSE(dCheck(1, 0), 4.0, tol);
  BOOST_CHECK_CLOSE(dCheck(1, 1), 5.0, tol);
  BOOST_CHECK_CLOSE(dCheck(1, 2), 6.0, tol);
  BOOST_CHECK_CLOSE(dCheck(2, 0), 7.0, tol);
  BOOST_CHECK_CLOSE(dCheck(2, 1), 8.0, tol);
  BOOST_CHECK_CLOSE(dCheck(2, 2), 3.0, tol);

  eig_decomp.factorize(sA, sEig, sR, sRinv);
  sCheck = sR * sEig * sRinv;

  BOOST_CHECK_CLOSE(sCheck(0, 0), 1.0, tol);
  BOOST_CHECK_CLOSE(sCheck(0, 1), 2.0, tol);
  BOOST_CHECK_CLOSE(sCheck(0, 2), 3.0, tol);
  BOOST_CHECK_CLOSE(sCheck(1, 0), 4.0, tol);
  BOOST_CHECK_CLOSE(sCheck(1, 1), 5.0, tol);
  BOOST_CHECK_CLOSE(sCheck(1, 2), 6.0, tol);
  BOOST_CHECK_CLOSE(sCheck(2, 0), 7.0, tol);
  BOOST_CHECK_CLOSE(sCheck(2, 1), 8.0, tol);
  BOOST_CHECK_CLOSE(sCheck(2, 2), 3.0, tol);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(LU_decomposition_utest)
{
  math::LUDecomposition<Real> lu_decomp;

  // lu_decomp.set_rank(rank);

  lu_decomp.factorize(dA, dL, dU);
  std::cout << "L = " << std::endl << dL << std::endl;
  std::cout << "U = " << std::endl << dU << std::endl;
  dCheck = dL * dU;

  const math::DenseConstVecView<Int> reordering = lu_decomp.reordering();

  std::cout << "A  = " << std::endl << dA << std::endl;

  std::cout << "L * U = " << std::endl << dCheck << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
