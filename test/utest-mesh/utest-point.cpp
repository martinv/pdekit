/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mesh_point_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <iostream>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "mesh/Point.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(point_constructor_test)
{
  /// Default constructor
  mesh::Point<_2D, Real> point2d;
  mesh::Point<_3D, Real> point3d;

  BOOST_CHECK_EQUAL(point2d.dim(), 2u);
  BOOST_CHECK_EQUAL(point3d.dim(), 3u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(point_indexing_test)
{
  // Default constructor
  mesh::Point<_2D, Real> point2d;
  mesh::Point<_3D, Real> point3d;

  point2d[X0] = 1.0;
  point2d[X1] = 2.0;

  BOOST_CHECK_EQUAL(point2d[X0], 1.0);
  BOOST_CHECK_EQUAL(point2d[X1], 2.0);

  point3d[X0] = 1.0;
  point3d[X1] = 2.0;
  point3d[X2] = 3.0;

  BOOST_CHECK_EQUAL(point3d[X0], 1.0);
  BOOST_CHECK_EQUAL(point3d[X1], 2.0);
  BOOST_CHECK_EQUAL(point3d[X2], 3.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(point_copy_assign_test)
{
  mesh::Point<_2D, Real> point2d, point2d_copy;
  mesh::Point<_3D, Real> point3d, point3d_copy;

  // Check for 2D points

  point2d[X0] = 1.0;
  point2d[X1] = 2.0;

  // Copy constructor
  point2d_copy = point2d;

  BOOST_CHECK_EQUAL(point2d_copy[X0], 1.0);
  BOOST_CHECK_EQUAL(point2d_copy[X1], 2.0);

  // Assignment operator
  point2d = point2d_copy;

  BOOST_CHECK_EQUAL(point2d[X0], 1.0);
  BOOST_CHECK_EQUAL(point2d[X1], 2.0);

  // Check for 3D points

  point3d[X0] = 1.0;
  point3d[X1] = 2.0;
  point3d[X2] = 3.0;

  point3d_copy = point3d;

  BOOST_CHECK_EQUAL(point3d_copy[X0], 1.0);
  BOOST_CHECK_EQUAL(point3d_copy[X1], 2.0);
  BOOST_CHECK_EQUAL(point3d_copy[X2], 3.0);

  point3d = point3d_copy;

  BOOST_CHECK_EQUAL(point3d[X0], 1.0);
  BOOST_CHECK_EQUAL(point3d[X1], 2.0);
  BOOST_CHECK_EQUAL(point3d[X2], 3.0);
}

// ----------------------------------------------------------------------------
