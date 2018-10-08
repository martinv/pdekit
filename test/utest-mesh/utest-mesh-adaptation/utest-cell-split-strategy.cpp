/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cell_split_strategy_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <ctime>
#include <iostream>

#include "mesh/CellTransform.hpp"
#include "mesh/adaptation/CellAdaptOp.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(utest_triag_uniform_refine)
{
  adapt::CellAdaptOp triag_adapt_op;
  triag_adapt_op.change_type(ElemShape::Triag, CellTransform::UNIFORM_REFINE);
  // triag_adapt_op.get().print();

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_parent_facets(), 3u);
  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_child_elems(), 4u);

  const PointSetTag expected_child_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  for (Uint c = 0; c < triag_adapt_op.get().nb_child_elems(); ++c)
  {
    const PointSetTag child_type(triag_adapt_op.get().child_elem_shape(c),
                                 expected_child_type.poly_order(),
                                 expected_child_type.ref_topology());
    BOOST_CHECK_EQUAL(child_type, expected_child_type);
  }

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_internal_child_facets(), 3u);

  const std::tuple<SUint, SUint> invalid_facet(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri0 = {
      std::tuple<SUint, SUint>(0, 0), invalid_facet, std::tuple<SUint, SUint>(2, 1)};

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri1 = {
      std::tuple<SUint, SUint>(0, 1), std::tuple<SUint, SUint>(1, 0), invalid_facet};

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri2 = {
      invalid_facet,
      std::tuple<SUint, SUint>(1, 1),
      std::tuple<SUint, SUint>(2, 0),
  };

  for (SUint f = 0; f < 3; ++f)
  {
    const std::tuple<SUint, SUint> ret0 = triag_adapt_op.get().containing_parent_facet_id(0u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret0), std::get<0>(ref_values_tri0[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret0), std::get<1>(ref_values_tri0[f]));

    const std::tuple<SUint, SUint> ret1 = triag_adapt_op.get().containing_parent_facet_id(1u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret1), std::get<0>(ref_values_tri1[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret1), std::get<1>(ref_values_tri1[f]));

    const std::tuple<SUint, SUint> ret2 = triag_adapt_op.get().containing_parent_facet_id(2u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret2), std::get<0>(ref_values_tri2[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret2), std::get<1>(ref_values_tri2[f]));

    // The fourth triangle is completely hidden inside the the parent
    // element and its 'child' facets are not incident to any facets of the
    // parent
    const std::tuple<SUint, SUint> ret3 = triag_adapt_op.get().containing_parent_facet_id(3u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret3), INVALID_LOC_ENTITY_ID);
    BOOST_CHECK_EQUAL(std::get<1>(ret3), INVALID_LOC_ENTITY_ID);
  }

  std::vector<math::DenseDMat<Real>> child_coords;

  PointSetTag child_tag(ElemShape::Triag, P5, PointSetID::Warpblend);
  triag_adapt_op.get().compute_child_coords(child_tag, child_coords);

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    std::cout << "Child coords (child " << c << ")" << std::endl;
    std::cout << child_coords[c] << std::endl << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(utest_triag_aniso_refine0)
{
  adapt::CellAdaptOp triag_adapt_op;
  triag_adapt_op.change_type(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_0);
  // triag_adapt_op.get().print();

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_parent_facets(), 3u);
  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_child_elems(), 2u);

  const PointSetTag expected_child_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  for (Uint c = 0; c < triag_adapt_op.get().nb_child_elems(); ++c)
  {
    const PointSetTag child_type(triag_adapt_op.get().child_elem_shape(c),
                                 expected_child_type.poly_order(),
                                 expected_child_type.ref_topology());
    BOOST_CHECK_EQUAL(child_type, expected_child_type);
  }

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_internal_child_facets(), 1u);

  const std::tuple<SUint, SUint> invalid_facet(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri0 = {
      std::tuple<SUint, SUint>(0, 0), invalid_facet, std::tuple<SUint, SUint>(2, 0)};

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri1 = {
      std::tuple<SUint, SUint>(0, 1), std::tuple<SUint, SUint>(1, 0), invalid_facet};

  for (SUint f = 0; f < 3; ++f)
  {
    const std::tuple<SUint, SUint> ret0 = triag_adapt_op.get().containing_parent_facet_id(0u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret0), std::get<0>(ref_values_tri0[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret0), std::get<1>(ref_values_tri0[f]));

    const std::tuple<SUint, SUint> ret1 = triag_adapt_op.get().containing_parent_facet_id(1u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret1), std::get<0>(ref_values_tri1[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret1), std::get<1>(ref_values_tri1[f]));
  }

  std::vector<math::DenseDMat<Real>> child_coords;

  PointSetTag child_tag(ElemShape::Triag, P5, PointSetID::Warpblend);
  triag_adapt_op.get().compute_child_coords(child_tag, child_coords);

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    std::cout << "Child coords (child " << c << ")" << std::endl;
    std::cout << child_coords[c] << std::endl << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(utest_triag_aniso_refine1)
{
  adapt::CellAdaptOp triag_adapt_op;
  triag_adapt_op.change_type(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_1);
  // triag_adapt_op.get().print();

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_parent_facets(), 3u);
  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_child_elems(), 2u);

  const PointSetTag expected_child_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  for (Uint c = 0; c < triag_adapt_op.get().nb_child_elems(); ++c)
  {
    const PointSetTag child_type(triag_adapt_op.get().child_elem_shape(c),
                                 expected_child_type.poly_order(),
                                 expected_child_type.ref_topology());
    BOOST_CHECK_EQUAL(child_type, expected_child_type);
  }

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_internal_child_facets(), 1u);

  const std::tuple<SUint, SUint> invalid_facet(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri0 = {
      std::tuple<SUint, SUint>(0, 0), std::tuple<SUint, SUint>(1, 0), invalid_facet};

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri1 = {
      invalid_facet, std::tuple<SUint, SUint>(1, 1), std::tuple<SUint, SUint>(2, 0)};

  for (SUint f = 0; f < 3; ++f)
  {
    const std::tuple<SUint, SUint> ret0 = triag_adapt_op.get().containing_parent_facet_id(0u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret0), std::get<0>(ref_values_tri0[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret0), std::get<1>(ref_values_tri0[f]));

    const std::tuple<SUint, SUint> ret1 = triag_adapt_op.get().containing_parent_facet_id(1u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret1), std::get<0>(ref_values_tri1[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret1), std::get<1>(ref_values_tri1[f]));
  }

  std::vector<math::DenseDMat<Real>> child_coords;

  PointSetTag child_tag(ElemShape::Triag, P5, PointSetID::Warpblend);
  triag_adapt_op.get().compute_child_coords(child_tag, child_coords);

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    std::cout << "Child coords (child " << c << ")" << std::endl;
    std::cout << child_coords[c] << std::endl << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(utest_triag_aniso_refine2)
{
  adapt::CellAdaptOp triag_adapt_op;
  triag_adapt_op.change_type(ElemShape::Triag, CellTransform::ANISO_REFINE_ORTHO_FACE_2);
  // triag_adapt_op.get().print();

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_parent_facets(), 3u);
  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_child_elems(), 2u);

  const PointSetTag expected_child_type = PointSetTag(ElemShape::Triag, P1, PointSetID::Equidist);

  for (Uint c = 0; c < triag_adapt_op.get().nb_child_elems(); ++c)
  {
    const PointSetTag child_type(triag_adapt_op.get().child_elem_shape(c),
                                 expected_child_type.poly_order(),
                                 expected_child_type.ref_topology());
    BOOST_CHECK_EQUAL(child_type, expected_child_type);
  }

  BOOST_CHECK_EQUAL(triag_adapt_op.get().nb_internal_child_facets(), 1u);

  const std::tuple<SUint, SUint> invalid_facet(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri0 = {
      std::tuple<SUint, SUint>(0, 0), invalid_facet, std::tuple<SUint, SUint>(2, 1)};

  const std::vector<std::tuple<SUint, SUint>> ref_values_tri1 = {
      invalid_facet, std::tuple<SUint, SUint>(1, 0), std::tuple<SUint, SUint>(2, 0)};

  for (SUint f = 0; f < 3; ++f)
  {
    const std::tuple<SUint, SUint> ret0 = triag_adapt_op.get().containing_parent_facet_id(0u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret0), std::get<0>(ref_values_tri0[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret0), std::get<1>(ref_values_tri0[f]));

    const std::tuple<SUint, SUint> ret1 = triag_adapt_op.get().containing_parent_facet_id(1u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret1), std::get<0>(ref_values_tri1[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret1), std::get<1>(ref_values_tri1[f]));
  }

  std::vector<math::DenseDMat<Real>> child_coords;

  PointSetTag child_tag(ElemShape::Triag, P5, PointSetID::Warpblend);
  triag_adapt_op.get().compute_child_coords(child_tag, child_coords);

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    std::cout << "Child coords (child " << c << ")" << std::endl;
    std::cout << child_coords[c] << std::endl << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(utest_quad_uniform_refine)
{
  adapt::CellAdaptOp quad_adapt_op;
  quad_adapt_op.change_type(ElemShape::Quad, CellTransform::UNIFORM_REFINE);

  BOOST_CHECK_EQUAL(quad_adapt_op.get().nb_parent_facets(), 4u);
  BOOST_CHECK_EQUAL(quad_adapt_op.get().nb_child_elems(), 4u);

  const PointSetTag expected_child_type = PointSetTag(ElemShape::Quad, P1, PointSetID::Equidist);

  for (Uint c = 0; c < quad_adapt_op.get().nb_child_elems(); ++c)
  {
    const PointSetTag child_type(quad_adapt_op.get().child_elem_shape(c),
                                 expected_child_type.poly_order(),
                                 expected_child_type.ref_topology());
    BOOST_CHECK_EQUAL(child_type, expected_child_type);
  }

  BOOST_CHECK_EQUAL(quad_adapt_op.get().nb_internal_child_facets(), 4u);

  const std::tuple<SUint, SUint> invalid_facet(INVALID_LOC_ENTITY_ID, INVALID_LOC_ENTITY_ID);

  const std::vector<std::tuple<SUint, SUint>> ref_values_quad0 = {
      std::tuple<SUint, SUint>(0, 0), invalid_facet, invalid_facet, std::tuple<SUint, SUint>(3, 1)};

  const std::vector<std::tuple<SUint, SUint>> ref_values_quad1 = {
      std::tuple<SUint, SUint>(0, 1), std::tuple<SUint, SUint>(1, 0), invalid_facet, invalid_facet};

  const std::vector<std::tuple<SUint, SUint>> ref_values_quad2 = {
      invalid_facet, std::tuple<SUint, SUint>(1, 1), std::tuple<SUint, SUint>(2, 0), invalid_facet};

  const std::vector<std::tuple<SUint, SUint>> ref_values_quad3 = {
      invalid_facet, invalid_facet, std::tuple<SUint, SUint>(2, 1), std::tuple<SUint, SUint>(3, 0)};

  for (SUint f = 0; f < 4; ++f)
  {
    const std::tuple<SUint, SUint> ret0 = quad_adapt_op.get().containing_parent_facet_id(0u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret0), std::get<0>(ref_values_quad0[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret0), std::get<1>(ref_values_quad0[f]));

    const std::tuple<SUint, SUint> ret1 = quad_adapt_op.get().containing_parent_facet_id(1u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret1), std::get<0>(ref_values_quad1[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret1), std::get<1>(ref_values_quad1[f]));

    const std::tuple<SUint, SUint> ret2 = quad_adapt_op.get().containing_parent_facet_id(2u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret2), std::get<0>(ref_values_quad2[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret2), std::get<1>(ref_values_quad2[f]));

    const std::tuple<SUint, SUint> ret3 = quad_adapt_op.get().containing_parent_facet_id(3u, f);
    BOOST_CHECK_EQUAL(std::get<0>(ret3), std::get<0>(ref_values_quad3[f]));
    BOOST_CHECK_EQUAL(std::get<1>(ret3), std::get<1>(ref_values_quad3[f]));
  }

  std::vector<math::DenseDMat<Real>> child_coords;

  PointSetTag child_tag(ElemShape::Quad, P5, PointSetID::Warpblend);
  quad_adapt_op.get().compute_child_coords(child_tag, child_coords);

  for (Uint c = 0; c < child_coords.size(); ++c)
  {
    std::cout << "Child coords (child " << c << ")" << std::endl;
    std::cout << child_coords[c] << std::endl << std::endl;
  }
}
