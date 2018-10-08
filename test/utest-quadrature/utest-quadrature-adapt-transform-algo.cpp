/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_adapt_transform_algo_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/Constants.hpp"
#include "common/PDEKit.hpp"
#include "mesh/point_set/QuadratureAdaptTransformAlgoFactory.hpp"

using namespace pdekit;
using namespace pdekit::mesh;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(line_quadrature_adapt_transform_utest)
{
  const Real tol = 1.e-13;

  QuadratureAdaptTransformAlgoFactory::instance_type &quad_trans_algo_factory =
      QuadratureAdaptTransformAlgoFactory::instance();

  const QuadratureAdaptTransformAlgoFactory::instance_type::const_product_base_ptr ta =
      quad_trans_algo_factory.create(mesh::CellTransform::UNIFORM_REFINE);

  math::DenseDMat<Real> qcoord_in(4, _1D);

  qcoord_in(0, X0) = -0.8;
  qcoord_in(1, X0) = -0.3;
  qcoord_in(2, X0) = 0.3;
  qcoord_in(3, X0) = 0.8;

  math::DenseDVec<Real> qwgt_in(4);

  qwgt_in[0] = 0.25;
  qwgt_in[1] = 0.75;
  qwgt_in[2] = 0.75;
  qwgt_in[3] = 0.25;

  math::DenseDMat<Real> qcoord_out(4, _1D);
  math::DenseDVec<Real> qwgt_out(4);

  // Test transformation for local id 0
  ta->compute_transformed_coords(qcoord_in, 0, qcoord_out);
  ta->compute_transformed_weights(qwgt_in, 0, qwgt_out);

  BOOST_CHECK_CLOSE(qcoord_out(0, X0), -0.9, tol);
  BOOST_CHECK_CLOSE(qcoord_out(1, X0), -0.65, tol);
  BOOST_CHECK_CLOSE(qcoord_out(2, X0), -0.35, tol);
  BOOST_CHECK_CLOSE(qcoord_out(3, X0), -0.1, tol);

  BOOST_CHECK_EQUAL(qwgt_in.size(), qwgt_out.size());

  for (Uint i = 0; i < qwgt_in.size(); ++i)
  {
    BOOST_CHECK_CLOSE(0.5 * qwgt_in[i], qwgt_out[i], tol);
  }

  // Test transformation for local id 1
  ta->compute_transformed_coords(qcoord_in, 1, qcoord_out);
  ta->compute_transformed_weights(qwgt_in, 1, qwgt_out);

  BOOST_CHECK_CLOSE(qcoord_out(0, X0), 0.1, tol);
  BOOST_CHECK_CLOSE(qcoord_out(1, X0), 0.35, tol);
  BOOST_CHECK_CLOSE(qcoord_out(2, X0), 0.65, tol);
  BOOST_CHECK_CLOSE(qcoord_out(3, X0), 0.9, tol);

  BOOST_CHECK_EQUAL(qwgt_in.size(), qwgt_out.size());

  for (Uint i = 0; i < qwgt_in.size(); ++i)
  {
    BOOST_CHECK_CLOSE(0.5 * qwgt_in[i], qwgt_out[i], tol);
  }

  /*
  std::cout << "Output quadrature coordinates: " << std::endl;
  std::cout << qcoord_out << std::endl;
  std::cout << "Output quadrature weights: " << std::endl;
  std::cout << qwgt_out << std::endl;
  */
}

// ----------------------------------------------------------------------------
