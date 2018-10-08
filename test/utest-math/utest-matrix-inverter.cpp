/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE matrix_inversion_utest
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "math/DenseDMat.hpp"
#include "math/DenseSMat.hpp"

using namespace pdekit;
using namespace pdekit::math;

struct MatrixInversionUtestFixture
{
  clock_t start, end;
  Real elapsed;
  static const Uint n_loops = 1e7;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MatrixInversion_TestSuite, MatrixInversionUtestFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_inverter_2x2)
{

  //  clock_t start, end;
  //  Real elapsed;
  //  const Uint n_loops = 1e7;

  std::cout << "-------------------------------------------------------------" << std::endl;
  std::cout << "\t2x2 MATRIX INVERSION" << std::endl;
  std::cout << "-------------------------------------------------------------" << std::endl;
  DenseSMat<Real, 2, 2> m2, m_inv2, unit2;

  m2(0, 0) = 1.0;
  m2(0, 1) = 2.0;
  m2(1, 0) = 3.0;
  m2(1, 1) = 4.0;

  start = clock();
  for (Uint loop = 0; loop < n_loops; ++loop)
  {
    m2.inv(m_inv2);
  }
  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "CPU time (2x2 inverter) = " << elapsed << " s" << std::endl;

  unit2 = m2 * m_inv2;

  const Real eps = 1.e-9;

  BOOST_CHECK_CLOSE(unit2(0, 0), 1.0, eps);
  BOOST_CHECK_CLOSE(unit2(1, 1), 1.0, eps);
  BOOST_CHECK_CLOSE(unit2(0, 1), 0.0, eps);
  BOOST_CHECK_CLOSE(unit2(1, 0), 0.0, eps);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_inverter_3x3)
{

  std::cout << "-------------------------------------------------------------" << std::endl;
  std::cout << "\t3x3 MATRIX INVERSION" << std::endl;
  std::cout << "-------------------------------------------------------------" << std::endl;

  DenseSMat<Real, 3, 3> m3, m_inv3, unit3;

  int sign = 1.0;

  for (Uint i = 0; i < 3; ++i)
  {
    for (Uint j = 0; j < 3; ++j)
    {
      sign *= -1.0;
      m3(i, j) = sign * (i * 3 + j + 1.0);
    }
  }
  m3(2, 2) = 4.0;

  start = clock();
  for (Uint loop = 0; loop < n_loops; ++loop)
  {
    m3.inv(m_inv3);
  }
  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "CPU time (3x3 inverter) = " << elapsed << " s" << std::endl;

  std::cout << "Check for 3x3 inversion:" << std::endl;
  unit3 = m3 * m_inv3;

  /*
  std::cout << "Reference solution (Matlab)" << std::endl;
  std::cout << "| 1.743589743589743   0.820512820512820   0.076923076923077 |"
  << std::endl; std::cout << "| 1.487179487179487   0.641025641025641
  0.153846153846154 |" << std::endl; std::cout << "| 0.076923076923077
  0.153846153846154   0.076923076923077 |" << std::endl
            << std::endl;
  */

  std::cout << "Value of determinant: det(M) = " << m3.det() << " (reference value = -39.0)"
            << std::endl
            << std::endl;

  const Real eps = 1.e-12;

  BOOST_CHECK_LE(std::abs(m3.det() + 39.0), eps);

  BOOST_CHECK_LE(std::abs(unit3(0, 0) - 1.0), eps);
  BOOST_CHECK_LE(std::abs(unit3(0, 1)), eps);
  BOOST_CHECK_LE(std::abs(unit3(0, 2)), eps);

  BOOST_CHECK_LE(std::abs(unit3(1, 0)), eps);
  BOOST_CHECK_LE(std::abs(unit3(1, 1) - 1.0), eps);
  BOOST_CHECK_LE(std::abs(unit3(1, 2)), eps);

  BOOST_CHECK_LE(std::abs(unit3(2, 0)), eps);
  BOOST_CHECK_LE(std::abs(unit3(2, 1)), eps);
  BOOST_CHECK_LE(std::abs(unit3(2, 2) - 1.0), eps);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_inverter_4x4)
{

  std::cout << "-------------------------------------------------------------" << std::endl;
  std::cout << "\t4x4 MATRIX INVERSION" << std::endl;
  std::cout << "-------------------------------------------------------------" << std::endl;

  DenseSMat<Real, 4, 4> m4, m_inv4, unit4;

  DenseDMat<Real> m4dyn, m_inv4dyn, unit4dyn;

  m4dyn.resize(4, 4);
  m_inv4dyn.resize(4, 4);
  unit4dyn.resize(4, 4);

  m4 = values_list(1.)(2.)(3.)(4.)(8.)(7.)(6.)(5.)(1.)(3.)(4.)(2.)(13.)(14.)(15.)(2.);

  m4dyn = m4;

  /// TIMING FOR STATIC MATRICES:
  start = clock();
  for (Uint loop = 0; loop < n_loops; ++loop)
  {
    m4.inv(m_inv4);
  }
  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "CPU time (4x4 static inverter) = " << elapsed << " s" << std::endl;

  /// TIMING FOR DYNAMIC MATRICES:
  start = clock();
  for (Uint loop = 0; loop < n_loops; ++loop)
  {
    m4dyn.inv(m_inv4dyn);
  }
  end     = clock();
  elapsed = ((Real)(end - start)) / CLOCKS_PER_SEC;

  std::cout.setf(std::ios::fixed);
  std::cout.precision(5);
  std::cout << "CPU time (4x4 dynamic inverter) = " << elapsed << " s" << std::endl;

  std::cout << "Check for 4x4 inversion:" << std::endl;
  unit4 = m4 * m_inv4;

  /*
  std::cout << "Reference solution (Matlab)" << std::endl;
  std::cout << "|  0.611111111111111  -0.174603174603175  -1.000000000000000 "
               "0.214285714285714 |" << std::endl;
  std::cout << "| -1.722222222222223   0.777777777777778   2.000000000000001 "
               "-0.500000000000000 |" << std::endl;
  std::cout << "|  1.055555555555556  -0.587301587301587  -1.000000000000001 "
               "0.357142857142857 |" << std::endl;
  std::cout << "|  0.166666666666667   0.095238095238095   0.000000000000000 "
               "-0.071428571428571 |" << std::endl
            << std::endl;
  */

  std::cout << "Value of determinant: det(M) = " << m4.det() << " (reference value = -126.0)"
            << std::endl
            << std::endl;

  const Real eps = 1.e-12;

  BOOST_CHECK_LE(std::abs(m4.det() + 126.0), eps);

  BOOST_CHECK_LE(std::abs(unit4(0, 0) - 1.0), eps);
  BOOST_CHECK_LE(std::abs(unit4(0, 1)), eps);
  BOOST_CHECK_LE(std::abs(unit4(0, 2)), eps);
  BOOST_CHECK_LE(std::abs(unit4(0, 3)), eps);

  BOOST_CHECK_LE(std::abs(unit4(1, 0)), eps);
  BOOST_CHECK_LE(std::abs(unit4(1, 1) - 1.0), eps);
  BOOST_CHECK_LE(std::abs(unit4(1, 2)), eps);
  BOOST_CHECK_LE(std::abs(unit4(1, 3)), eps);

  BOOST_CHECK_LE(std::abs(unit4(2, 0)), eps);
  BOOST_CHECK_LE(std::abs(unit4(2, 1)), eps);
  BOOST_CHECK_LE(std::abs(unit4(2, 2) - 1.0), eps);
  BOOST_CHECK_LE(std::abs(unit4(2, 3)), eps);

  BOOST_CHECK_LE(std::abs(unit4(3, 0)), eps);
  BOOST_CHECK_LE(std::abs(unit4(3, 1)), eps);
  BOOST_CHECK_LE(std::abs(unit4(3, 2)), eps);
  BOOST_CHECK_LE(std::abs(unit4(3, 3) - 1.0), eps);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(matrix_inverter_6x6)
{

  DenseSMat<Real, 6, 6> a;
  DenseSMat<Real, 6, 6> a_inverse;

  DenseSMat<Real, 6, 6> unit6;

  a(0, 0) = 1.0;
  a(1, 1) = 2.0;
  a(2, 2) = 3.0;
  a(3, 3) = 4.0;
  a(4, 4) = 5.0;
  a(5, 2) = 200;
  a(5, 5) = -1.;

  a(0, 4) = 20.0;
  a(1, 2) = -3.0;
  a(2, 0) = 7.0;

  a.inv(a_inverse);

  unit6 = a * a_inverse;

  std::cout << "Product a*a_inv = " << std::endl;
  std::cout << unit6 << std::endl;

  const Real eps = 5.e-15;

  for (Uint i = 0; i < 6; ++i)
  {
    BOOST_CHECK_LE(std::abs(unit6(i, i) - 1.0), eps);
  }

  for (Uint i = 0; i < 6; ++i)
    for (Uint j = 0; j < 6; ++j)
    {
      if (i != j)
      {
        BOOST_CHECK_LE(std::abs(unit6(i, j)), eps);
      }
    }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
