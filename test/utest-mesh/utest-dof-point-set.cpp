/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE dof_point_set_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "common/PDEKit.hpp"
#include "math/DenseDVec.hpp"
#include "math/DenseSMat.hpp"
#include "math/MathConstants.hpp"
#include "mesh/point_set/QuadratureTransformUtils.hpp"
#include "mesh/point_set/StdPointSet.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::common;

// ----------------------------------------------------------------------------

Real sum_weights(const math::DenseDVec<Real> &weights)
{
  Real sum = 0.0;
  for (Uint i = 0; i < weights.size(); ++i)
  {
    sum += weights[i];
  }
  return sum;
}

// ----------------------------------------------------------------------------

void check_weight_sum(const ElemShape elem_shape, const PointSetID pt_set_id, const Uint deg_low,
                      const Uint deg_high, const Real ref_sum)
{
  StdPointSet qd;

  for (Uint poly_order = deg_low; poly_order <= deg_high; ++poly_order)
  {
    qd.change_type(elem_shape, poly_order, pt_set_id);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, ref_sum, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

void check_triag_point_set_rotation_permutation(const math::DenseDMat<Real> &quad_coords,
                                                const std::vector<Uint> &rot_permutation_vec,
                                                const Real eps_tol = 1.e-14)
{
  math::DenseDMat<Real> bar_coords(quad_coords.rows(), quad_coords.cols());
  math::DenseDMat<Real> equ_coords(quad_coords.rows(), quad_coords.cols());
  math::DenseDMat<Real> rotated_equ_coords(quad_coords.rows(), quad_coords.cols());

  // Vertices of reference triangle
  math::DenseSVec<Real, 2> V0, V1, V2;
  V0[X0] = -1.0;
  V0[X1] = -1.0;

  V1[X0] = 1.0;
  V1[X1] = -1.0;

  V2[X0] = -1.0;
  V2[X1] = 1.0;

  mesh::detail::compute_triag_barycentric_coords(V0, V1, V2, quad_coords, bar_coords);

  // Vertices of an equilateral triangle WHOSE CENTER IS AT (0, 0)
  // (If the center is at (0,0), it's easier to rotate the triangle
  math::DenseSVec<Real, 2> V0_equ, V1_equ, V2_equ;
  V0_equ[X0] = -1.0;
  V0_equ[X1] = -1. / 3. * std::sqrt(3.);

  V1_equ[X0] = 1.0;
  V1_equ[X1] = -1. / 3. * std::sqrt(3.);

  V2_equ[X0] = 0.0;
  V2_equ[X1] = 2. / 3. * std::sqrt(3.);

  const Uint nb_qd_pts = quad_coords.rows();

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    equ_coords(qd_pt, X0) = V0_equ[X0] * bar_coords(qd_pt, 0) + V1_equ[X0] * bar_coords(qd_pt, 1) +
                            V2_equ[X0] * bar_coords(qd_pt, 2);
    equ_coords(qd_pt, X1) = V0_equ[X1] * bar_coords(qd_pt, 0) + V1_equ[X1] * bar_coords(qd_pt, 1) +
                            V2_equ[X1] * bar_coords(qd_pt, 2);
  }

  // Transformation matrix for rotation by 120 degrees
  math::DenseSMat<Real, 2, 2> rot_matrix;
  const Real angle = 120. * math::pi / 180.;
  rot_matrix(0, 0) = std::cos(angle);
  rot_matrix(0, 1) = -std::sin(angle);
  rot_matrix(1, 0) = std::sin(angle);
  rot_matrix(1, 1) = std::cos(angle);

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    rotated_equ_coords(qd_pt, X0) =
        rot_matrix(0, 0) * equ_coords(qd_pt, X0) + rot_matrix(0, 1) * equ_coords(qd_pt, X1);
    rotated_equ_coords(qd_pt, X1) =
        rot_matrix(1, 0) * equ_coords(qd_pt, X0) + rot_matrix(1, 1) * equ_coords(qd_pt, X1);
  }

  std::cout << "Original coordinates: " << std::endl << equ_coords << std::endl;
  std::cout << "Rotated coordinates: " << std::endl << rotated_equ_coords << std::endl;
  std::cout << "Permutation vector:" << std::endl;
  for (Uint i = 0; i < rot_permutation_vec.size(); ++i)
  {
    std::cout << rot_permutation_vec[i] << " ";
  }
  std::cout << std::endl;

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    BOOST_CHECK_LE(
        std::abs(rotated_equ_coords(rot_permutation_vec[qd_pt], X0) - equ_coords(qd_pt, X0)),
        eps_tol);
    BOOST_CHECK_LE(
        std::abs(rotated_equ_coords(rot_permutation_vec[qd_pt], X1) - equ_coords(qd_pt, X1)),
        eps_tol);
  }
}

// ----------------------------------------------------------------------------

void check_triag_point_set_flip_permutation(const math::DenseDMat<Real> &quad_coords,
                                            const std::vector<Uint> &flip_permutation_vec,
                                            const Real eps_tol = 1.e-14)
{
  math::DenseDMat<Real> bar_coords(quad_coords.rows(), quad_coords.cols());
  math::DenseDMat<Real> equ_coords(quad_coords.rows(), quad_coords.cols());
  math::DenseDMat<Real> flipped_equ_coords(quad_coords.rows(), quad_coords.cols());

  // Vertices of reference triangle
  math::DenseSVec<Real, 2> V0, V1, V2;
  V0[X0] = -1.0;
  V0[X1] = -1.0;

  V1[X0] = 1.0;
  V1[X1] = -1.0;

  V2[X0] = -1.0;
  V2[X1] = 1.0;

  mesh::detail::compute_triag_barycentric_coords(V0, V1, V2, quad_coords, bar_coords);

  // Vertices of a triangle which is symmetric with respect to the y-axis
  math::DenseSVec<Real, 2> V0_equ, V1_equ, V2_equ;
  V0_equ[X0] = -1.0;
  V0_equ[X1] = -1.5;

  V1_equ[X0] = 1.0;
  V1_equ[X1] = -1.5;

  V2_equ[X0] = 0.0;
  V2_equ[X1] = 1.0;

  const Uint nb_qd_pts = quad_coords.rows();

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    equ_coords(qd_pt, X0) = V0_equ[X0] * bar_coords(qd_pt, 0) + V1_equ[X0] * bar_coords(qd_pt, 1) +
                            V2_equ[X0] * bar_coords(qd_pt, 2);
    equ_coords(qd_pt, X1) = V0_equ[X1] * bar_coords(qd_pt, 0) + V1_equ[X1] * bar_coords(qd_pt, 1) +
                            V2_equ[X1] * bar_coords(qd_pt, 2);
  }

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    flipped_equ_coords(qd_pt, X0) = -1. * equ_coords(qd_pt, X0);
    // y-coordinate does not change
    flipped_equ_coords(qd_pt, X1) = equ_coords(qd_pt, X1);
  }

  // std::cout << "Original coordinates: " << std::endl << equ_coords <<
  // std::endl; std::cout << "Flipped coordinates: " << std::endl <<
  // flipped_equ_coords << std::endl;

  for (Uint qd_pt = 0; qd_pt < nb_qd_pts; ++qd_pt)
  {
    BOOST_CHECK_LE(
        std::abs(flipped_equ_coords(flip_permutation_vec[qd_pt], X0) - equ_coords(qd_pt, X0)),
        eps_tol);
    BOOST_CHECK_LE(
        std::abs(flipped_equ_coords(flip_permutation_vec[qd_pt], X1) - equ_coords(qd_pt, X1)),
        eps_tol);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_line_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Line, PointSetID::Equidist, P1, P15, 2.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_line_warpblend_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Line, PointSetID::Warpblend, P1, P15, 2.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_triag_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Triag, PointSetID::Equidist, P1, P15, 2.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_triag_warpblend_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Triag, PointSetID::Warpblend, P1, P15, 2.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_quad_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Quad, PointSetID::Equidist, P1, P10, 4.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_quad_warpblend_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Quad, PointSetID::Warpblend, P1, P10, 4.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_tetra_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Tetra, PointSetID::Equidist, P1, P9, 4.0 / 3.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_tetra_warpblend_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Tetra, PointSetID::Warpblend, P1, P9, 4.0 / 3.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_hexa_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Hexa, PointSetID::Equidist, P1, P6, 8.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(std_pyramid_equidist_pt_set_weights_utest)
{
  check_weight_sum(ElemShape::Pyramid, PointSetID::Equidist, P1, P6, 8.0 / 3.0);
}

// ----------------------------------------------------------------------------
