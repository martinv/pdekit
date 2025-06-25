/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE range_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "common/Range1D.hpp"

using namespace pdekit;
using namespace pdekit::common;

BOOST_AUTO_TEST_CASE(range_utest)
{
  Range1D<Uint> range1(2, 5);

  Range1D<Uint> range2 = range1;
  range2 += 3;

  BOOST_CHECK_EQUAL(range1.lbound(), 2u);
  BOOST_CHECK_EQUAL(range1.ubound(), 5u);

  BOOST_CHECK_EQUAL(range2.lbound(), 5u);
  BOOST_CHECK_EQUAL(range2.ubound(), 8u);

  BOOST_CHECK_EQUAL(range1.size(), 4u);
  BOOST_CHECK_EQUAL(range2.size(), 4u);

  BOOST_CHECK_EQUAL(range1.in_range(4), true);
  BOOST_CHECK_EQUAL(range1.in_range(9), false);

  std::cout << range1 << std::endl;
  std::cout << range2 << std::endl;

  const Range1D<Uint> range3 = {2, 3};

  BOOST_CHECK_EQUAL(range3.lbound(), 2u);
  BOOST_CHECK_EQUAL(range3.ubound(), 3u);

  const Range1D<Uint> range4 = {3, 3};
  const Range1D<Uint> range5 = {4, 3};

  BOOST_CHECK_EQUAL(range4.empty(), false);
  BOOST_CHECK_EQUAL(range5.empty(), true);
}
