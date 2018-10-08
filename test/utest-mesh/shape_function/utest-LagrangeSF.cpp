/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Lagrange_sf_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "common/PDEKit.hpp"
#include "fixture/ReferenceShapeFunctionsTetra.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/shape_function/ShapeFunction.hpp"
#include "mesh/std_region/StdRegion.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::common;

// ----------------------------------------------------------------------------

struct LagrangeSFUtest_Fixture
{
  /// common setup for each test case
  LagrangeSFUtest_Fixture();

  /// common tear-down for each test case
  ~LagrangeSFUtest_Fixture();

  /// Check that the sum of elements in each matrix row
  /// is equal to given reference value with specified tolerance
  void check_row_sums(const math::DenseDMat<Real> &mat, const Real ref, const Real tol) const;

  /// Compute polynomial values in points
  void compute_poly_values_tri(const Uint P, const math::DenseDMat<Real> &points,
                               math::DenseDVec<Real> &values) const;

  /// Compute polynomial derivatives in points
  void compute_poly_derivatives_tri(const Uint P, const math::DenseDMat<Real> &points,
                                    math::DenseDMat<Real> &values) const;

  /// Compute polynomial values in points
  void compute_poly_values_tet(const Uint P, const math::DenseDMat<Real> &points,
                               math::DenseDVec<Real> &values) const;

  /// Compute polynomial derivatives in points
  void compute_poly_derivatives_tet(const Uint P, const math::DenseDMat<Real> &points,
                                    math::DenseDMat<Real> &values) const;

  // P1 polynomials on triangle
  const Real poly_triag_p1(const Real xi0, const Real xi1) const;
  const Real deriv_xi0_poly_triag_p1(const Real xi0, const Real xi1) const;
  const Real deriv_xi1_poly_triag_p1(const Real xi0, const Real xi1) const;

  // P2 polynomials on triangle
  const Real poly_triag_p2(const Real xi0, const Real xi1) const;
  const Real deriv_xi0_poly_triag_p2(const Real xi0, const Real xi1) const;
  const Real deriv_xi1_poly_triag_p2(const Real xi0, const Real xi1) const;

  // P3 polynomials on triangle
  const Real poly_triag_p3(const Real xi0, const Real xi1) const;
  const Real deriv_xi0_poly_triag_p3(const Real xi0, const Real xi1) const;
  const Real deriv_xi1_poly_triag_p3(const Real xi0, const Real xi1) const;

  // P1 polynomials on tetrahedron
  const Real poly_tet_p1(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi0_poly_tet_p1(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi1_poly_tet_p1(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi2_poly_tet_p1(const Real xi0, const Real xi1, const Real xi2) const;

  // P2 polynomials on tetrahedron
  const Real poly_tet_p2(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi0_poly_tet_p2(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi1_poly_tet_p2(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi2_poly_tet_p2(const Real xi0, const Real xi1, const Real xi2) const;

  // P3 polynomials on tetrahedron
  const Real poly_tet_p3(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi0_poly_tet_p3(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi1_poly_tet_p3(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi2_poly_tet_p3(const Real xi0, const Real xi1, const Real xi2) const;

  // P4 polynomials on tetrahedron
  const Real poly_tet_p4(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi0_poly_tet_p4(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi1_poly_tet_p4(const Real xi0, const Real xi1, const Real xi2) const;
  const Real deriv_xi2_poly_tet_p4(const Real xi0, const Real xi1, const Real xi2) const;
};

LagrangeSFUtest_Fixture::LagrangeSFUtest_Fixture()
{
}

LagrangeSFUtest_Fixture::~LagrangeSFUtest_Fixture()
{
}

void LagrangeSFUtest_Fixture::check_row_sums(const math::DenseDMat<Real> &mat, const Real ref,
                                             const Real tol) const
{
  for (Uint r = 0; r < mat.rows(); ++r)
  {
    Real sum = 0.0;
    for (Uint c = 0; c < mat.cols(); ++c)
    {
      sum += mat(r, c);
    }

    const Real diff = sum - ref;

    BOOST_CHECK_LE(std::abs(diff), tol);
  }
}

void LagrangeSFUtest_Fixture::compute_poly_values_tri(const Uint P,
                                                      const math::DenseDMat<Real> &points,
                                                      math::DenseDVec<Real> &values) const
{
  values.resize(points.rows());

  switch (P)
  {
    case 1:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_triag_p1(points(i, X0), points(i, X1));
      }
      break;
    }
    case 2:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_triag_p2(points(i, X0), points(i, X1));
      }
      break;
    }
    case 3:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_triag_p3(points(i, X0), points(i, X1));
      }
      break;
    }
  }
}

void LagrangeSFUtest_Fixture::compute_poly_derivatives_tri(const Uint P,
                                                           const math::DenseDMat<Real> &points,
                                                           math::DenseDMat<Real> &values) const
{
  values.resize(points.rows(), _2D);

  switch (P)
  {
    case 1:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_triag_p1(points(i, X0), points(i, X1));
        values(i, X1) = deriv_xi1_poly_triag_p1(points(i, X0), points(i, X1));
      }
      break;
    }
    case 2:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_triag_p2(points(i, X0), points(i, X1));
        values(i, X1) = deriv_xi1_poly_triag_p2(points(i, X0), points(i, X1));
      }
      break;
    }
    case 3:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_triag_p3(points(i, X0), points(i, X1));
        values(i, X1) = deriv_xi1_poly_triag_p3(points(i, X0), points(i, X1));
      }
      break;
    }
  }
}

void LagrangeSFUtest_Fixture::compute_poly_values_tet(const Uint P,
                                                      const math::DenseDMat<Real> &points,
                                                      math::DenseDVec<Real> &values) const
{
  values.resize(points.rows());

  switch (P)
  {
    case 1:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_tet_p1(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 2:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_tet_p2(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 3:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_tet_p3(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 4:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values[i] = poly_tet_p4(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
  }
}

void LagrangeSFUtest_Fixture::compute_poly_derivatives_tet(const Uint P,
                                                           const math::DenseDMat<Real> &points,
                                                           math::DenseDMat<Real> &values) const
{
  values.resize(points.rows(), _3D);

  switch (P)
  {
    case 1:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_tet_p1(points(i, X0), points(i, X1), points(i, X2));
        values(i, X1) = deriv_xi1_poly_tet_p1(points(i, X0), points(i, X1), points(i, X2));
        values(i, X2) = deriv_xi2_poly_tet_p1(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 2:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_tet_p2(points(i, X0), points(i, X1), points(i, X2));
        values(i, X1) = deriv_xi1_poly_tet_p2(points(i, X0), points(i, X1), points(i, X2));
        values(i, X2) = deriv_xi2_poly_tet_p2(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 3:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_tet_p3(points(i, X0), points(i, X1), points(i, X2));
        values(i, X1) = deriv_xi1_poly_tet_p3(points(i, X0), points(i, X1), points(i, X2));
        values(i, X2) = deriv_xi2_poly_tet_p3(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
    case 4:
    {
      for (Uint i = 0; i < points.rows(); ++i)
      {
        values(i, X0) = deriv_xi0_poly_tet_p4(points(i, X0), points(i, X1), points(i, X2));
        values(i, X1) = deriv_xi1_poly_tet_p4(points(i, X0), points(i, X1), points(i, X2));
        values(i, X2) = deriv_xi2_poly_tet_p4(points(i, X0), points(i, X1), points(i, X2));
      }
      break;
    }
  }
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_triag_p1(const Real xi0, const Real xi1) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_triag_p1(const Real xi0, const Real xi1) const
{
  return 4.2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_triag_p1(const Real xi0, const Real xi1) const
{
  return -9.3;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_triag_p2(const Real xi0, const Real xi1) const
{

  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 5.3 * xi0 * xi1 + 2.2 * xi0 * xi0 - 8.4 * xi1 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_triag_p2(const Real xi0, const Real xi1) const
{
  return 4.2 + 5.3 * xi1 + 2.0 * 2.2 * xi0;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_triag_p2(const Real xi0, const Real xi1) const
{
  return -9.3 + 5.3 * xi0 - 2.0 * 8.4 * xi1;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_triag_p3(const Real xi0, const Real xi1) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 5.3 * xi0 * xi1 + 2.2 * xi0 * xi0 - 8.4 * xi1 * xi1 +
         2.1 * xi0 * xi0 * xi0 + 5.6 * xi0 * xi0 * xi1 - 2.2 * xi0 * xi1 * xi1 +
         0.456 * xi1 * xi1 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_triag_p3(const Real xi0, const Real xi1) const
{
  return 4.2 + 5.3 * xi1 + 2.0 * 2.2 * xi0 + 3.0 * 2.1 * xi0 * xi0 + 2.0 * 5.6 * xi0 * xi1 -
         2.2 * xi1 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_triag_p3(const Real xi0, const Real xi1) const
{
  return -9.3 + 5.3 * xi0 - 2.0 * 8.4 * xi1 + 5.6 * xi0 * xi0 - 2.0 * 2.2 * xi0 * xi1 +
         3.0 * 0.456 * xi1 * xi1;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_tet_p1(const Real xi0, const Real xi1,
                                                const Real xi2) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 0.185 * xi2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_tet_p1(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 4.2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_tet_p1(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return -9.3;
}

const Real LagrangeSFUtest_Fixture::deriv_xi2_poly_tet_p1(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 0.185;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_tet_p2(const Real xi0, const Real xi1,
                                                const Real xi2) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 0.185 * xi2 + 1.84 * xi0 * xi1 + 5.893 * xi0 * xi2 +
         9.23 * xi1 * xi2 + 1.23 * xi0 * xi0 - 9.3 * xi1 * xi1 + 0.32 * xi2 * xi2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_tet_p2(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 4.2 + 1.84 * xi1 + 5.893 * xi2 + 2.0 * 1.23 * xi0;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_tet_p2(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return -9.3 + 1.84 * xi0 + 9.23 * xi2 - 2.0 * 9.3 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi2_poly_tet_p2(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 0.185 + 5.893 * xi0 + 9.23 * xi1 + 2.0 * 0.32 * xi2;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_tet_p3(const Real xi0, const Real xi1,
                                                const Real xi2) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 0.185 * xi2 + 1.84 * xi0 * xi1 + 5.893 * xi0 * xi2 +
         9.23 * xi1 * xi2 + 1.23 * xi0 * xi0 - 9.3 * xi1 * xi1 + 0.32 * xi2 * xi2 +
         2.12 * xi0 * xi1 * xi2 + 4.43 * xi0 * xi0 * xi1 + 3.36 * xi0 * xi1 * xi1 -
         2.34154 * xi0 * xi0 * xi2 - 5.35632 * xi0 * xi2 * xi2 - 8.313 * xi1 * xi1 * xi2 +
         3.845 * xi1 * xi2 * xi2 + 1.1234 * xi0 * xi0 * xi0 - 2.234 * xi1 * xi1 * xi1 +
         6.245 * xi2 * xi2 * xi2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_tet_p3(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 4.2 + 1.84 * xi1 + 5.893 * xi2 + 2.0 * 1.23 * xi0 + 2.12 * xi1 * xi2 +
         2.0 * 4.43 * xi0 * xi1 + 3.36 * xi1 * xi1 - 2.0 * 2.34154 * xi0 * xi2 -
         5.35632 * xi2 * xi2 + 3.0 * 1.1234 * xi0 * xi0;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_tet_p3(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return -9.3 + 1.84 * xi0 + 9.23 * xi2 - 2.0 * 9.3 * xi1 + 2.12 * xi0 * xi2 + 4.43 * xi0 * xi0 +
         2.0 * 3.36 * xi0 * xi1 - 2.0 * 8.313 * xi1 * xi2 + 3.845 * xi2 * xi2 -
         3.0 * 2.234 * xi1 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi2_poly_tet_p3(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 0.185 + 5.893 * xi0 + 9.23 * xi1 + 2.0 * 0.32 * xi2 + 2.12 * xi0 * xi1 -
         2.34154 * xi0 * xi0 - 2.0 * 5.35632 * xi0 * xi2 - 8.313 * xi1 * xi1 +
         2.0 * 3.845 * xi1 * xi2 + 3.0 * 6.245 * xi2 * xi2;
}

// ----------------------------------------------------------------------------

const Real LagrangeSFUtest_Fixture::poly_tet_p4(const Real xi0, const Real xi1,
                                                const Real xi2) const
{
  return 3.2 + 4.2 * xi0 - 9.3 * xi1 + 0.185 * xi2 + 1.84 * xi0 * xi1 + 5.893 * xi0 * xi2 +
         9.23 * xi1 * xi2 + 1.23 * xi0 * xi0 - 9.3 * xi1 * xi1 + 0.32 * xi2 * xi2 +
         2.12 * xi0 * xi1 * xi2 + 4.43 * xi0 * xi0 * xi1 + 3.36 * xi0 * xi1 * xi1 -
         2.34154 * xi0 * xi0 * xi2 - 5.35632 * xi0 * xi2 * xi2 - 8.313 * xi1 * xi1 * xi2 +
         3.845 * xi1 * xi2 * xi2 + 1.1234 * xi0 * xi0 * xi0 - 2.234 * xi1 * xi1 * xi1 +
         6.245 * xi2 * xi2 * xi2 + 1.253 * xi0 * xi1 * xi1 * xi1 - 1.276 * xi0 * xi0 * xi0 * xi1 +
         3.273 * xi0 * xi2 * xi2 * xi2 + 5.676 * xi0 * xi0 * xi0 * xi2 +
         7.243 * xi1 * xi2 * xi2 * xi2 - 3.206 * xi1 * xi1 * xi1 * xi2 -
         4.41 * xi0 * xi0 * xi1 * xi2 + 6.542 * xi0 * xi1 * xi1 * xi2 +
         0.031 * xi0 * xi1 * xi2 * xi2 + 0.31 * xi0 * xi0 * xi2 * xi2 -
         6.41 * xi0 * xi0 * xi1 * xi1 + 0.041 * xi1 * xi1 * xi2 * xi2 +
         1.534 * xi0 * xi0 * xi0 * xi0 + 8.1004 * xi1 * xi1 * xi1 * xi1 -
         0.00341 * xi2 * xi2 * xi2 * xi2;
}

const Real LagrangeSFUtest_Fixture::deriv_xi0_poly_tet_p4(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 4.2 + 1.84 * xi1 + 5.893 * xi2 + 2.0 * 1.23 * xi0 + 2.12 * xi1 * xi2 +
         2.0 * 4.43 * xi0 * xi1 + 3.36 * xi1 * xi1 - 2.0 * 2.34154 * xi0 * xi2 -
         5.35632 * xi2 * xi2 + 3.0 * 1.1234 * xi0 * xi0 + 1.253 * xi1 * xi1 * xi1 -
         3.0 * 1.276 * xi0 * xi0 * xi1 + 3.273 * xi2 * xi2 * xi2 + 3.0 * 5.676 * xi0 * xi0 * xi2 -
         2.0 * 4.41 * xi0 * xi1 * xi2 + 6.542 * xi1 * xi1 * xi2 + 0.031 * xi1 * xi2 * xi2 +
         2.0 * 0.31 * xi0 * xi2 * xi2 - 2.0 * 6.41 * xi0 * xi1 * xi1 +
         4.0 * 1.534 * xi0 * xi0 * xi0;
}

const Real LagrangeSFUtest_Fixture::deriv_xi1_poly_tet_p4(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return -9.3 + 1.84 * xi0 + 9.23 * xi2 - 2.0 * 9.3 * xi1 + 2.12 * xi0 * xi2 + 4.43 * xi0 * xi0 +
         2.0 * 3.36 * xi0 * xi1 - 2.0 * 8.313 * xi1 * xi2 + 3.845 * xi2 * xi2 -
         3.0 * 2.234 * xi1 * xi1 + 3.0 * 1.253 * xi0 * xi1 * xi1 - 1.276 * xi0 * xi0 * xi0 +
         7.243 * xi2 * xi2 * xi2 - 3.0 * 3.206 * xi1 * xi1 * xi2 - 4.41 * xi0 * xi0 * xi2 +
         2.0 * 6.542 * xi0 * xi1 * xi2 + 0.031 * xi0 * xi2 * xi2 - 2.0 * 6.41 * xi0 * xi0 * xi1 +
         2.0 * 0.041 * xi1 * xi2 * xi2 + 4.0 * 8.1004 * xi1 * xi1 * xi1;
}

const Real LagrangeSFUtest_Fixture::deriv_xi2_poly_tet_p4(const Real xi0, const Real xi1,
                                                          const Real xi2) const
{
  return 0.185 + 5.893 * xi0 + 9.23 * xi1 + 2.0 * 0.32 * xi2 + 2.12 * xi0 * xi1 -
         2.34154 * xi0 * xi0 - 2.0 * 5.35632 * xi0 * xi2 - 8.313 * xi1 * xi1 +
         2.0 * 3.845 * xi1 * xi2 + 3.0 * 6.245 * xi2 * xi2 + 3.0 * 3.273 * xi0 * xi2 * xi2 +
         5.676 * xi0 * xi0 * xi0 + 3.0 * 7.243 * xi1 * xi2 * xi2 - 3.206 * xi1 * xi1 * xi1 -
         4.41 * xi0 * xi0 * xi1 + 6.542 * xi0 * xi1 * xi1 + 2.0 * 0.031 * xi0 * xi1 * xi2 +
         2.0 * 0.31 * xi0 * xi0 * xi2 + 2.0 * 0.041 * xi1 * xi1 * xi2 -
         4.0 * 0.00341 * xi2 * xi2 * xi2;
}

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(LagrangeSF_TestSuite, LagrangeSFUtest_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Lagrange_sf_utest)
{
  sf::detail::ShapeFunctionInstance sf_instance;
  const PointSetTag std_reg_tag(ElemShape::Triag, P1, PointSetID::Equidist);
  const sf::SFTag sf_tag(ElemShape::Triag, SFunc::Lagrange, P1, ModalBasis::Modal);
  sf::detail::ShapeFunctionInstance::construct(std::make_tuple(std_reg_tag, sf_tag), sf_instance);

  sf::ShapeFunction sf;
  sf.change_type(std_reg_tag, sf_tag);

  math::DenseDMat<Real> ref_coords(3, 2);
  ref_coords(0, 0) = -1.0;
  ref_coords(0, 1) = -1.0;

  ref_coords(1, 0) = 0.0;
  ref_coords(1, 1) = -0.3;

  ref_coords(2, 0) = 1.0;
  ref_coords(2, 1) = -1.0;

  math::DenseDMat<Real> sf_values;
  std::vector<math::DenseDMat<Real>> sf_derivatives;

  sf.get().compute_ref_values(ref_coords, sf_values);
  sf.get().compute_ref_derivatives(ref_coords, sf_derivatives);

  std::cout << "Points at which shape function values are computed:" << std::endl;
  std::cout << ref_coords << std::endl;

  std::cout << "Shape function values:" << std::endl;
  std::cout << sf_values << std::endl;

  const std::string CoordNames[3] = {"X", "Y", "Z"};

  for (Uint d = 0; d < sf_derivatives.size(); ++d)
  {
    std::cout << "Derivatives of shape functions with respect to " << CoordNames[d] << ":"
              << std::endl;
    std::cout << sf_derivatives[d] << std::endl;
  }

  // sf.change_type(sf::ShapeFunctionInstance::tag_type(Triag,Lagrange,Nodal,P5,Dubiner));

  // sf::NewShapeFunction newsf;

  // ==========================================================================
  // Verify that the shape function values and derivatives are identical when
  // computed using the primal basis and using hardcoded values from utest
  // fixture
  // ==========================================================================

  const PointSetTag std_reg_tag_tet(ElemShape::Tetra, P3, PointSetID::Equidist);
  const sf::SFTag sf_tag_tet(ElemShape::Tetra, SFunc::Lagrange, P3, ModalBasis::Modal);

  sf.change_type(std_reg_tag_tet, sf_tag_tet);
  mesh::StdPointSet quad;
  quad.change_type(ElemShape::Tetra, P5, PointSetID::Gauss);

  math::DenseDMat<Real> computed_Lagrange_values;
  math::DenseDMat<Real> reference_Lagrange_values;
  std::vector<math::DenseDMat<Real>> computed_Lagrange_derivatives(_3D);
  std::vector<math::DenseDMat<Real>> reference_Lagrange_derivatives(_3D);
  math::DenseDMat<Real> value_differences;
  math::DenseDVec<Real> reference_point;
  math::DenseDVec<Real> sf_values_in_one_point;
  math::DenseDMat<Real> sf_derivatives_in_one_point;

  // Compute Lagrange shape function values and their derivatives
  sf.get().compute_ref_values(quad.get().coordinates(), computed_Lagrange_values);
  sf.get().compute_ref_derivatives(quad.get().coordinates(), computed_Lagrange_derivatives);

  const Uint nb_sf    = computed_Lagrange_values.cols();
  const Uint nb_qd_pt = computed_Lagrange_values.rows();

  reference_point.resize(_3D);
  sf_values_in_one_point.resize(nb_sf);
  sf_derivatives_in_one_point.resize(nb_sf, _3D);

  reference_Lagrange_values.resize(nb_qd_pt, nb_sf);
  for (Uint d = 0; d < _3D; ++d)
  {
    reference_Lagrange_derivatives[d].resize(nb_qd_pt, nb_sf);
  }

  for (Uint q = 0; q < nb_qd_pt; ++q)
  {
    reference_point[XI0] = quad.get().coordinates()(q, XI0);
    reference_point[XI1] = quad.get().coordinates()(q, XI1);
    reference_point[XI2] = quad.get().coordinates()(q, XI2);

    utest_fixture::p3_tetra_Lagrange_sf::eval(reference_point, sf_values_in_one_point);
    reference_Lagrange_values.insert_row(q, sf_values_in_one_point);

    for (Uint d = 0; d < _3D; ++d)
    {
      utest_fixture::p3_tetra_Lagrange_sf::eval_deriv(reference_point, sf_derivatives_in_one_point);
      reference_Lagrange_derivatives[d].insert_row(q, sf_derivatives_in_one_point.const_col(d));
    }
  }

  std::cout << "Computed Lagrange sf values: " << std::endl;
  std::cout << computed_Lagrange_values << std::endl;
  std::cout << "reference_Lagrange_values:" << std::endl;
  std::cout << reference_Lagrange_values << std::endl;

  value_differences.resize(nb_qd_pt, nb_sf);

  value_differences = computed_Lagrange_values - reference_Lagrange_values;

  Real sum;

  for (Uint r = 0; r < value_differences.rows(); ++r)
  {
    sum = 0.0;

    for (Uint c = 0; c < value_differences.cols(); ++c)
    {
      sum += computed_Lagrange_values(r, c);
      if (std::abs(value_differences(r, c)) < 1.e-15)
      {
        value_differences(r, c) = 0.0;
      }
    }

    // Check that the sum of Lagrange shape functions in one point is one
    BOOST_CHECK_CLOSE(sum, 1.0,
                      1.e-13); // HERE THE TOLERANCE IS IN PERCENT!
  }

  std::cout << "Difference (computed - reference) shape function values:" << std::endl;
  std::cout << value_differences << std::endl;

  for (Uint d = 0; d < _3D; ++d)
  {
    value_differences = computed_Lagrange_derivatives[d] - reference_Lagrange_derivatives[d];

    for (Uint r = 0; r < value_differences.rows(); ++r)
    {
      sum = 0.0;
      for (Uint c = 0; c < value_differences.cols(); ++c)
      {
        sum += computed_Lagrange_derivatives[d](r, c);
        if (std::abs(value_differences(r, c)) < 1.e-14)
        {
          value_differences(r, c) = 0.0;
        }
      }

      // Check that the sum of Lagrange sf. derivatives in one point is
      // zero
      BOOST_CHECK_LE(std::abs(sum), 1.e-13);
    }

    std::cout << "Difference (computed - reference) shape function derivatives "
                 "(with respect to var. "
              << d << " :" << std::endl;
    std::cout << value_differences << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Lagrange_sf_triag_utest)
{
  // Generate points at which the shape functions should be evaluated

  const Real x0 = -1.0;
  const Real y0 = -1.0;

  const Real x1 = 1.0;
  const Real y1 = -1.0;

  const Real x2 = -1.0;
  const Real y2 = 1.0;

  const Uint N     = 14;
  const Real delta = 1.0 / (N - 1);

  const Uint n_values = N * (N + 1) / 2;

  math::DenseDMat<Real> pt_coords;
  math::DenseDMat<Real> Vandermonde;
  std::vector<math::DenseDMat<Real>> Vandermonde_deriv;

  math::DenseDVec<Real> u_at_dofs;

  math::DenseDVec<Real> u_at_pts_check;
  math::DenseDMat<Real> du_at_pts_check;

  pt_coords.resize(n_values, 2);

  u_at_pts_check.resize(n_values);
  du_at_pts_check.resize(n_values, _2D);

  Uint pt_idx = 0;

  for (Uint i = 0; i < N; ++i)
  {
    const Real l0 = i * delta;

    for (Uint j = 0; j < N - i; ++j)
    {
      const Real l1 = j * delta;
      const Real l2 = 1.0 - l0 - l1;

      pt_coords(pt_idx, X0) = l0 * x0 + l1 * x1 + l2 * x2;
      pt_coords(pt_idx, X1) = l0 * y0 + l1 * y1 + l2 * y2;

      pt_idx++;
    }
  }

  sf::ShapeFunction sf_tri;
  mesh::StdRegion std_region_tri;

  for (Uint p_order = 1; p_order <= P3; ++p_order)
  {
    const mesh::PointSetTag tri_std_reg_tag(ElemShape::Triag, p_order, PointSetID::Equidist);
    std_region_tri.change_type(tri_std_reg_tag);

    // std::cout << "P = " << p_order << std::endl;
    const sf::SFTag tri_sf_tag(ElemShape::Triag, SFunc::Lagrange, p_order, ModalBasis::Modal);

    sf_tri.change_type(tri_std_reg_tag, tri_sf_tag);

    sf_tri.get().compute_ref_values(pt_coords, Vandermonde);
    sf_tri.get().compute_ref_derivatives(pt_coords, Vandermonde_deriv);

    check_row_sums(Vandermonde, 1.0, 1.e-14);
    check_row_sums(Vandermonde_deriv[0], 0.0, 1.e-14);
    check_row_sums(Vandermonde_deriv[1], 0.0, 1.e-14);

    u_at_dofs.resize(std_region_tri.get().nb_nodes());

    math::DenseDMat<Real> const &dof_coords = std_region_tri.get().coordinates();

    compute_poly_values_tri(p_order, dof_coords, u_at_dofs);

    compute_poly_values_tri(p_order, pt_coords, u_at_pts_check);
    compute_poly_derivatives_tri(p_order, pt_coords, du_at_pts_check);

    /*
    std::cout << "Values at dofs: " << std::endl;
    for (Uint i = 0; i < Vandermonde.cols(); ++i)
    {
      std::cout << u_at_dofs[i] << std::endl;
    }

    std::cout << "Values at check points: " << std::endl;
    for (Uint i = 0; i < u_at_pts_check.size(); ++i)
    {
      std::cout << u_at_pts_check[i] << std::endl;
    }
    */

    for (Uint i = 0; i < n_values; ++i)
    {
      Real u_h     = 0.0;
      Real du_h_dx = 0.0;
      Real du_h_dy = 0.0;

      for (Uint idof = 0; idof < Vandermonde.cols(); ++idof)
      {
        u_h += Vandermonde(i, idof) * u_at_dofs[idof];
        du_h_dx += Vandermonde_deriv[X0](i, idof) * u_at_dofs[idof];
        du_h_dy += Vandermonde_deriv[X1](i, idof) * u_at_dofs[idof];
      }

      // std::cout << "u = " << u_at_pts_check[i] << ", u_h = " << u_h <<
      // std::endl;

      const Real delta_u     = std::abs(u_h - u_at_pts_check[i]);
      const Real delta_du_dx = std::abs(du_h_dx - du_at_pts_check(i, X0));
      const Real delta_du_dy = std::abs(du_h_dy - du_at_pts_check(i, X1));
      BOOST_CHECK_LE(delta_u, 1.e-13);
      BOOST_CHECK_LE(delta_du_dx, 1.e-13);
      BOOST_CHECK_LE(delta_du_dy, 1.e-13);
    }
  }

  // std::cout << sf_values << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Lagrange_sf_tetra_utest)
{
  // Generate points at which the shape functions should be evaluated

  const Real x0 = -1.0;
  const Real y0 = -1.0;
  const Real z0 = -1.0;

  const Real x1 = 1.0;
  const Real y1 = -1.0;
  const Real z1 = -1.0;

  const Real x2 = -1.0;
  const Real y2 = 1.0;
  const Real z2 = -1.0;

  const Real x3 = -1.0;
  const Real y3 = -1.0;
  const Real z3 = 1.0;

  const Uint N     = 14;
  const Real delta = 1.0 / (N - 1);

  const Uint n_values = N * (N + 1) * (N + 2) / 6;

  math::DenseDMat<Real> pt_coords;
  math::DenseDMat<Real> Vandermonde;
  std::vector<math::DenseDMat<Real>> Vandermonde_deriv;

  math::DenseDVec<Real> u_at_dofs;

  math::DenseDVec<Real> u_at_pts_check;
  math::DenseDMat<Real> du_at_pts_check;

  pt_coords.resize(n_values, 3);

  u_at_pts_check.resize(n_values);
  du_at_pts_check.resize(n_values, _3D);

  Uint pt_idx = 0;

  for (Uint i = 0; i < N; ++i)
  {
    const Real l0 = i * delta;

    for (Uint j = 0; j < N - i; ++j)
    {
      const Real l1 = j * delta;

      for (Uint k = 0; k < N - i - j; ++k)
      {
        const Real l2 = k * delta;
        const Real l3 = 1.0 - l0 - l1 - l2;

        pt_coords(pt_idx, X0) = l0 * x0 + l1 * x1 + l2 * x2 + l3 * x3;
        pt_coords(pt_idx, X1) = l0 * y0 + l1 * y1 + l2 * y2 + l3 * y3;
        pt_coords(pt_idx, X2) = l0 * z0 + l1 * z1 + l2 * z2 + l3 * z3;
        pt_idx++;
      }
    }
  }

  sf::ShapeFunction sf_tet;
  mesh::StdRegion std_region_tet;

  for (Uint p_order = 1; p_order <= P4; ++p_order)
  {
    const mesh::PointSetTag tet_std_reg_tag(ElemShape::Tetra, p_order, PointSetID::Equidist);
    std_region_tet.change_type(tet_std_reg_tag);

    // std::cout << "P = " << p_order << std::endl;
    const sf::SFTag tet_sf_tag(ElemShape::Tetra, SFunc::Lagrange, p_order, ModalBasis::Modal);

    sf_tet.change_type(tet_std_reg_tag, tet_sf_tag);

    sf_tet.get().compute_ref_values(pt_coords, Vandermonde);
    sf_tet.get().compute_ref_derivatives(pt_coords, Vandermonde_deriv);

    check_row_sums(Vandermonde, 1.0, 1.e-14);
    check_row_sums(Vandermonde_deriv[0], 0.0, 1.e-14);
    check_row_sums(Vandermonde_deriv[1], 0.0, 1.e-14);
    check_row_sums(Vandermonde_deriv[2], 0.0, 1.e-14);

    u_at_dofs.resize(std_region_tet.get().nb_nodes());

    math::DenseDMat<Real> const &dof_coords = std_region_tet.get().coordinates();

    compute_poly_values_tet(p_order, dof_coords, u_at_dofs);

    compute_poly_values_tet(p_order, pt_coords, u_at_pts_check);
    compute_poly_derivatives_tet(p_order, pt_coords, du_at_pts_check);

    /*
    std::cout << "Values at dofs: " << std::endl;
    for (Uint i = 0; i < Vandermonde.cols(); ++i)
    {
      std::cout << u_at_dofs[i] << std::endl;
    }

    std::cout << "Values at check points: " << std::endl;
    for (Uint i = 0; i < u_at_pts_check.size(); ++i)
    {
      std::cout << u_at_pts_check[i] << std::endl;
    }
    */

    for (Uint i = 0; i < n_values; ++i)
    {
      Real u_h     = 0.0;
      Real du_h_dx = 0.0;
      Real du_h_dy = 0.0;
      Real du_h_dz = 0.0;

      for (Uint idof = 0; idof < Vandermonde.cols(); ++idof)
      {
        u_h += Vandermonde(i, idof) * u_at_dofs[idof];
        du_h_dx += Vandermonde_deriv[X0](i, idof) * u_at_dofs[idof];
        du_h_dy += Vandermonde_deriv[X1](i, idof) * u_at_dofs[idof];
        du_h_dz += Vandermonde_deriv[X2](i, idof) * u_at_dofs[idof];
      }

      // std::cout << "u = " << u_at_pts_check[i] << ", u_h = " << u_h <<
      // std::endl;

      const Real delta_u     = std::abs(u_h - u_at_pts_check[i]);
      const Real delta_du_dx = std::abs(du_h_dx - du_at_pts_check(i, X0));
      const Real delta_du_dy = std::abs(du_h_dy - du_at_pts_check(i, X1));
      const Real delta_du_dz = std::abs(du_h_dz - du_at_pts_check(i, X2));
      BOOST_CHECK_LE(delta_u, 1.e-12);
      BOOST_CHECK_LE(delta_du_dx, 1.e-12);
      BOOST_CHECK_LE(delta_du_dy, 1.e-12);
      BOOST_CHECK_LE(delta_du_dz, 1.e-12);
    }
  }

  // std::cout << sf_values << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
