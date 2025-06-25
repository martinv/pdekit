/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_permutation_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/point_set/QuadraturePermutation.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(line_quadrature_permutation_utest)
{

  QuadraturePermutation quad_permutation;

  mesh::PointSetTag qtag(ElemShape::Line, P1, PointSetID::Gauss);
  mesh::EntityRealignCode perm_code;
  perm_code.add_flip();

  quad_permutation.change_type(qtag, 0, perm_code);

  BOOST_CHECK_EQUAL(quad_permutation.get().type_id().as_string(), "Line-P1-Gauss");
  BOOST_CHECK_EQUAL(quad_permutation.get().local_id(), 0u);

  BOOST_CHECK_EQUAL(quad_permutation.get().size(), 2u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(triag_quadrature_permutation_utest)
{
  QuadraturePermutation quad_permutation;

  mesh::PointSetTag qtag(ElemShape::Triag, P2, PointSetID::Gauss);
  mesh::EntityRealignCode perm_code;
  perm_code.add_flip();
  perm_code.add_rotation();

  quad_permutation.change_type(qtag, 1, perm_code);

  std::cout << "{" << qtag.as_string() << "} [";
  for (Uint i = 0; i < quad_permutation.get().size(); ++i)
  {
    std::cout << " " << quad_permutation.get().vertex(i);
  }
  std::cout << " ]" << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(triag_face_quadrature_permutation_utest)
{
  QuadraturePermutation quad_permutation;

  mesh::PointSetTag qtag(ElemShape::Triag, P1, PointSetID::FaceGauss);
  mesh::EntityRealignCode perm_code;
  perm_code.add_flip();

  quad_permutation.change_type(qtag, 1, perm_code);

  BOOST_CHECK_EQUAL(quad_permutation.get().type_id().as_string(), "Triag-P1-FaceGauss");
  BOOST_CHECK_EQUAL(quad_permutation.get().local_id(), 1u);
  BOOST_CHECK_EQUAL(quad_permutation.get().size(), 2u);
  BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 1u);
  BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 0u);

  // Quadratures for order 2 and 3
  for (Uint p = 2; p <= 3; ++p)
  {
    qtag                = mesh::PointSetTag(ElemShape::Triag, p, PointSetID::FaceGauss);
    const Uint local_id = p % 3; // Vary the local id a bit
    quad_permutation.change_type(qtag, local_id, perm_code);
    BOOST_CHECK_EQUAL(quad_permutation.get().size(), 2u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 1u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 0u);
  }

  // Quadratures for order 4 and 5
  for (Uint p = 4; p <= 5; ++p)
  {
    qtag                = mesh::PointSetTag(ElemShape::Triag, p, PointSetID::FaceGauss);
    const Uint local_id = p % 3; // Vary the local id a bit
    quad_permutation.change_type(qtag, local_id, perm_code);
    BOOST_CHECK_EQUAL(quad_permutation.get().size(), 3u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 2u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 1u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(2), 0u);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(quad_face_quadrature_permutation_utest)
{
  QuadraturePermutation quad_permutation;

  mesh::PointSetTag qtag(ElemShape::Quad, P1, PointSetID::FaceGauss);
  mesh::EntityRealignCode perm_code;
  perm_code.add_flip();

  quad_permutation.change_type(qtag, 1, perm_code);

  BOOST_CHECK_EQUAL(quad_permutation.get().type_id().as_string(), "Quad-P1-FaceGauss");
  BOOST_CHECK_EQUAL(quad_permutation.get().local_id(), 1u);
  BOOST_CHECK_EQUAL(quad_permutation.get().size(), 2u);
  BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 1u);
  BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 0u);

  // Quadratures for order 2 and 3
  // We will test for 1 flip and 2 rotations - should be the
  // same as a single flip (1D edge rotated twice is the same
  // edge we had at the beginning)
  perm_code.set_nb_flips(1);
  perm_code.set_nb_rotations(2);

  for (Uint p = 2; p <= 3; ++p)
  {
    qtag                = mesh::PointSetTag(ElemShape::Quad, p, PointSetID::FaceGauss);
    const Uint local_id = p % 4; // Vary the local id a bit
    quad_permutation.change_type(qtag, local_id, perm_code);
    BOOST_CHECK_EQUAL(quad_permutation.get().size(), 2u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 1u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 0u);
  }

  // Quadratures for order 4
  for (Uint p = 4; p <= 4; ++p)
  {
    qtag                = mesh::PointSetTag(ElemShape::Quad, p, PointSetID::FaceGauss);
    const Uint local_id = p % 4; // Vary the local id a bit
    quad_permutation.change_type(qtag, local_id, perm_code);
    BOOST_CHECK_EQUAL(quad_permutation.get().size(), 3u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(0), 2u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(1), 1u);
    BOOST_CHECK_EQUAL(quad_permutation.get().vertex(2), 0u);
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tetra_face_quadrature_permutation_utest)
{
  QuadraturePermutation quad_permutation;

  mesh::PointSetTag qtag(ElemShape::Tetra, P2, PointSetID::FaceGauss);
  mesh::EntityRealignCode perm_code;
  perm_code.add_flip();
  perm_code.add_rotation();

  quad_permutation.change_type(qtag, 1, perm_code);

  std::cout << "{" << qtag.as_string() << "} [";
  for (Uint i = 0; i < quad_permutation.get().size(); ++i)
  {
    std::cout << " " << quad_permutation.get().vertex(i);
  }
  std::cout << " ]" << std::endl;
}

// ----------------------------------------------------------------------------
