/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE array_view_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

#include "common/ArrayShape.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(array_view_construction)
{
  // Check default constructor
  common::ArrayShape<2, Uint> shape_2d_a;
  BOOST_CHECK_EQUAL(shape_2d_a.nb_dims(), 2u);

  // Check constructor for multiple arguments
  common::ArrayShape<5, Uint> shape_5d_a(2, 3, 4, 5, 6);
  BOOST_CHECK_EQUAL(shape_5d_a.nb_dims(), 5u);
  BOOST_CHECK_EQUAL(shape_5d_a.size(0), 2u);
  BOOST_CHECK_EQUAL(shape_5d_a.size(1), 3u);
  BOOST_CHECK_EQUAL(shape_5d_a.size(2), 4u);
  BOOST_CHECK_EQUAL(shape_5d_a.size(3), 5u);
  BOOST_CHECK_EQUAL(shape_5d_a.size(4), 6u);

  // Check copy constructor
  common::ArrayShape<5, Uint> shape_5d_b(shape_5d_a);
  BOOST_CHECK_EQUAL(shape_5d_b.nb_dims(), 5u);
  BOOST_CHECK_EQUAL(shape_5d_b.size(0), 2u);
  BOOST_CHECK_EQUAL(shape_5d_b.size(1), 3u);
  BOOST_CHECK_EQUAL(shape_5d_b.size(2), 4u);
  BOOST_CHECK_EQUAL(shape_5d_b.size(3), 5u);
  BOOST_CHECK_EQUAL(shape_5d_b.size(4), 6u);

  // Check assignment operator
  common::ArrayShape<5, Uint> shape_5d_c;
  shape_5d_c = shape_5d_b;

  BOOST_CHECK_EQUAL(shape_5d_c.nb_dims(), 5u);
  BOOST_CHECK_EQUAL(shape_5d_c.size(0), 2u);
  BOOST_CHECK_EQUAL(shape_5d_c.size(1), 3u);
  BOOST_CHECK_EQUAL(shape_5d_c.size(2), 4u);
  BOOST_CHECK_EQUAL(shape_5d_c.size(3), 5u);
  BOOST_CHECK_EQUAL(shape_5d_c.size(4), 6u);
}

// ----------------------------------------------------------------------------
