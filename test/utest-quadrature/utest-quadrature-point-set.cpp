/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_point_set_test
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

void check_triag_quadrature_rotation_permutation(const math::DenseDMat<Real> &quad_coords,
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

void check_triag_quadrature_flip_permutation(const math::DenseDMat<Real> &quad_coords,
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

BOOST_AUTO_TEST_CASE(Gauss_line_quadrature_weights_utest)
{
  StdPointSet qd;
  qd.change_type(ElemShape::Line, P1, PointSetID::Gauss);

  std::cout << "Coordinates of quadrature Line-P1-Gauss" << std::endl;
  std::cout << qd.get().coordinates() << std::endl;
  std::cout << "Weights of quadrature Line-P1-Gauss" << std::endl;
  std::cout << qd.get().weights() << std::endl;

  for (Uint poly_order = 1; poly_order <= 15; ++poly_order)
  {
    qd.change_type(ElemShape::Line, poly_order, PointSetID::Gauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, 2.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(GaussLobatto_line_quadrature_weights_utest)
{
  StdPointSet qd;
  qd.change_type(ElemShape::Line, P1, PointSetID::GaussLobatto);

  std::cout << "Coordinates of quadrature Line-P1-Gauss" << std::endl;
  std::cout << qd.get().coordinates() << std::endl;
  std::cout << "Weights of quadrature Line-P1-GaussLobatto" << std::endl;
  std::cout << qd.get().weights() << std::endl;

  for (Uint poly_order = 1; poly_order <= 15; ++poly_order)
  {
    qd.change_type(ElemShape::Line, poly_order, PointSetID::GaussLobatto);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, 2.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_triag_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 15; ++poly_order)
  {
    qd.change_type(ElemShape::Triag, poly_order, PointSetID::Gauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, 2.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_triag_face_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 5; ++poly_order)
  {
    qd.change_type(ElemShape::Triag, poly_order, PointSetID::FaceGauss);

    math::DenseDVec<Real> const &weights_face_0 = qd.get().weights(0);
    const Real sum_w0                           = sum_weights(weights_face_0);
    BOOST_CHECK_CLOSE(sum_w0, 2.0, 1.e-10);

    math::DenseDVec<Real> const &weights_face_1 = qd.get().weights(1);
    const Real sum_w1                           = sum_weights(weights_face_1);
    BOOST_CHECK_CLOSE(sum_w1, std::sqrt(2.) * 2.0, 1.e-10);

    math::DenseDVec<Real> const &weights_face_2 = qd.get().weights(2);
    const Real sum_w2                           = sum_weights(weights_face_2);
    BOOST_CHECK_CLOSE(sum_w2, 2.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_quad_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 8; ++poly_order)
  {
    qd.change_type(ElemShape::Quad, poly_order, PointSetID::Gauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, 4.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_quad_face_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 4; ++poly_order)
  {
    qd.change_type(ElemShape::Quad, poly_order, PointSetID::FaceGauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_CLOSE(sum_w, 2.0, 1.e-10);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_tetra_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 10; ++poly_order)
  {
    qd.change_type(ElemShape::Tetra, poly_order, PointSetID::Gauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_LE(std::abs(sum_w - 4.0 / 3.0), 3.e-8);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_tetra_face_quadrature_weights_utest)
{
  StdPointSet qd;

  const Real ref_sum_w_face012 = 2.0;
  const Real ref_sum_w_face3   = std::sqrt(3) * 2.0;

  for (Uint poly_order = 1; poly_order <= 5; ++poly_order)
  {
    qd.change_type(ElemShape::Tetra, poly_order, PointSetID::FaceGauss);

    math::DenseDVec<Real> const &weights_face_0 = qd.get().weights(0);
    const Real sum_w0                           = sum_weights(weights_face_0);
    BOOST_CHECK_CLOSE(sum_w0, ref_sum_w_face012, 1.e-10);

    math::DenseDVec<Real> const &weights_face_1 = qd.get().weights(1);
    const Real sum_w1                           = sum_weights(weights_face_1);
    BOOST_CHECK_CLOSE(sum_w1, ref_sum_w_face012, 1.e-10);

    math::DenseDVec<Real> const &weights_face_2 = qd.get().weights(2);
    const Real sum_w2                           = sum_weights(weights_face_2);
    BOOST_CHECK_CLOSE(sum_w2, ref_sum_w_face012, 1.e-10);

    math::DenseDVec<Real> const &weights_face_3 = qd.get().weights(3);
    const Real sum_w3                           = sum_weights(weights_face_3);
    BOOST_CHECK_CLOSE(sum_w3, ref_sum_w_face3, 1.e-10);
  }

  qd.change_type(ElemShape::Tetra, P5, PointSetID::FaceGauss);

  std::ofstream outfile;
  outfile.open("tetra_face_quad_pts.3D");
  outfile << "x y z face_idx" << std::endl;
  outfile.setf(std::ios::fixed);
  outfile.precision(15);

  for (Uint f = 0; f < qd.get().nb_local_entities(); ++f)
  {
    math::DenseDMat<Real> const &coords = qd.get().coordinates(f);

    for (Uint q = 0; q < coords.rows(); ++q)
    {
      outfile << std::setw(20) << coords(q, XI0) << " " << std::setw(20) << coords(q, XI1) << " "
              << std::setw(20) << coords(q, XI2) << " " << f << std::endl;
    }
  }

  outfile.close();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Gauss_hexa_quadrature_weights_utest)
{
  StdPointSet qd;

  for (Uint poly_order = 1; poly_order <= 7; ++poly_order)
  {
    qd.change_type(ElemShape::Hexa, poly_order, PointSetID::Gauss);
    math::DenseDVec<Real> const &weights = qd.get().weights();

    Real sum_w = sum_weights(weights);
    BOOST_CHECK_LE(std::abs(sum_w - 8.0), 3.e-8);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(compute_barycentric_coords)
{
  StdPointSet qd;

  qd.change_type(ElemShape::Triag, P5, PointSetID::Gauss);

  const math::DenseDMat<Real> &coords = qd.get().coordinates();

  math::DenseDMat<Real> bar_coords, new_cart_coords;

  // We suppose that the vertices of the reference triangle are
  // (-1,-1), (1,-1), (-1,1)
  math::DenseSVec<Real, 2> V0, V1, V2;
  V0[X0] = -1.;
  V0[X1] = -1.;

  V1[X0] = 1.;
  V1[X1] = -1.;

  V2[X0] = -1.;
  V2[X1] = 1.;

  mesh::detail::compute_triag_barycentric_coords(V0, V1, V2, coords, bar_coords);

  // std::cout << "Original coordinates = " << std::endl << coords <<
  // std::endl; std::cout << "Barycentric coordinates = " << std::endl <<
  // bar_coords << std::endl;

  new_cart_coords.resize(coords.rows(), coords.cols());

  // Try to rescale the barycentric coordinates to fit into
  // a triangle with vertices
  // (0,-2), (3, -2), (0, 1)

  for (Uint qd_pt = 0; qd_pt < new_cart_coords.rows(); ++qd_pt)
  {

    const Real x0 =
        0.0 * bar_coords(qd_pt, 0) + 3.0 * bar_coords(qd_pt, 1) + 0.0 * bar_coords(qd_pt, 2);
    const Real x1 =
        -2.0 * bar_coords(qd_pt, 0) - 2.0 * bar_coords(qd_pt, 1) + 1.0 * bar_coords(qd_pt, 2);
    new_cart_coords(qd_pt, X0) = x0;
    new_cart_coords(qd_pt, X1) = x1;
  }

  // std::cout << "New Cartesian coordinates: " << std::endl <<
  // new_cart_coords
  // << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(compute_triag_quadrature_rotation_permutation)
{
  StdPointSet qd;
  std::vector<Uint> rot_permutation_vec;

  for (Uint p = 1; p <= P15; ++p)
  {
    qd.change_type(ElemShape::Triag, p, PointSetID::Gauss);
    mesh::detail::canonical_triag_quadrature_rotation_permutation(qd.get().coordinates(),
                                                                  rot_permutation_vec);

    const Real tol_eps = p < P9 ? 1.e-14 : 1.e-11;
    check_triag_quadrature_rotation_permutation(qd.get().coordinates(), rot_permutation_vec,
                                                tol_eps);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(compute_triag_quadrature_flip_permutation)
{
  StdPointSet qd;
  std::vector<Uint> flip_permutation_vec;

  for (Uint p = 1; p <= P15; ++p)
  {
    qd.change_type(ElemShape::Triag, p, PointSetID::Gauss);
    mesh::detail::canonical_triag_quadrature_flip_permutation(qd.get().coordinates(),
                                                              flip_permutation_vec);
    const Real tol_eps = p < P9 ? 1.e-14 : 1.e-11;
    check_triag_quadrature_flip_permutation(qd.get().coordinates(), flip_permutation_vec, tol_eps);
  }
}

// ----------------------------------------------------------------------------
