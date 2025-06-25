/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE matrix_cond_number_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/DenseDMat.hpp"
#include "math/DenseSMat.hpp"
#include "math/unary_ops/MatrixCondNumber.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ============================================================================

BOOST_AUTO_TEST_CASE(matrix_condition_number_utest)
{
  math::DenseDMat<Real> A(3, 3);

  A = math::values_list(1.)(2.)(3.)(4.)(5.)(6.)(7.)(8.)(3.);

  const Real cn3_1norm   = math::cond_number_1_norm(A);
  const Real cn3_infnorm = math::cond_number_inf_norm(A);

  BOOST_CHECK_CLOSE(cn3_1norm, 55.0,
                    1.e-13); // Reference value provided by Matlab
  BOOST_CHECK_CLOSE(cn3_infnorm, 54.0,
                    1.e-13); // Reference value provided by Matlab

  std::cout << "3 x 3 matrix:" << std::endl;
  std::cout << "  Condition number (1-norm) = " << cn3_1norm << std::endl;
  std::cout << "  Condition number (infinity norm) = " << cn3_infnorm << std::endl;

  A.resize(4, 4);
  A = math::values_list(-1.)(2.)(3.)(4.)(5.)(-6.)(7.)(8.)(9.)(10.)(-11.)(12.)(13.)(14.)(15.)(-16.);

  const Real cn4_1norm   = math::cond_number_1_norm(A);
  const Real cn4_infnorm = math::cond_number_inf_norm(A);

  BOOST_CHECK_CLOSE(cn4_1norm, 17.231985940246044,
                    1.e-13); // Reference value provided by Matlab
  // BOOST_CHECK_CLOSE(cn4_infnorm, 14.168717047451668, 1.e-13); // Reference
  // by Matlab

  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "4 x 4 matrix:" << std::endl;
  std::cout << "  Condition number (1-norm) = " << cn4_1norm << std::endl;
  std::cout << "  Condition number (infinity norm) = " << cn4_infnorm << std::endl;

  A.resize(5, 5);
  A = math::values_list(-1.)(2.)(3.)(4.)(5.)(6.)(-7.)(8.)(9.)(10.)(11.)(12.)(-13.)(14.)(15.)(16.)(
      17.)(18.)(-19.)(20.)(21.)(22.)(23.)(24.)(-25.);

  const Real cn5_1norm   = math::cond_number_1_norm(A);
  const Real cn5_infnorm = math::cond_number_inf_norm(A);

  BOOST_CHECK_CLOSE(cn5_1norm, 30.990183456710657,
                    1.e-13); // Reference value provided by Matlab
  BOOST_CHECK_CLOSE(cn5_infnorm, 27.852429996781463,
                    1.e-13); // Reference value provided by Matlab

  std::cout << "A = " << std::endl << A << std::endl;
  std::cout << "5 x 5 matrix:" << std::endl;
  std::cout << "  Condition number (1-norm) = " << cn5_1norm << std::endl;
  std::cout << "  Condition number (infinity norm) = " << cn5_infnorm << std::endl;
}

// ============================================================================
