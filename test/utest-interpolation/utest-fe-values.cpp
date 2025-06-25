/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE fe_values_utest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>

#include "common/PDEKit.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/point_set/StdPointSetFactory.hpp"

using namespace pdekit;
using namespace pdekit::mesh;
using namespace pdekit::interpolation;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(fe_values_basic_test)
{
  PointSetTag std_region_tag(ElemShape::Triag, P1, PointSetID::Equidist);
  sf::SFTag type(ElemShape::Triag, SFunc::Lagrange, P1, ModalBasis::Modal);

  FEValues fe(std_region_tag, type);

  mesh::StdPointSet quadrature;
  quadrature.change_type(ElemShape::Triag, P2, PointSetID::Gauss);

  fe.fill_Vandermonde(quadrature.get().coordinates(), quadrature.get().weights());

  std::cout << "Vandermonde matrix:" << std::endl;
  std::cout << fe.Vandermonde() << std::endl;

  std::cout << "Vandermonde matrix of partial derivatives:" << std::endl;
  std::cout << fe.deriv_Vandermonde(X) << std::endl;
  std::cout << fe.deriv_Vandermonde(Y) << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(fe_values_filter_quad)
{
  PointSetTag std_region_tag(ElemShape::Quad, P3, PointSetID::Equidist);
  sf::SFTag sf_tag(ElemShape::Quad, SFunc::Lagrange, P3, ModalBasis::Modal);

  FEValues fe(std_region_tag, sf_tag);

  mesh::StdPointSet quad;
  quad.change_type(ElemShape::Quad, P3, PointSetID::Gauss);

  fe.fill_Vandermonde(quad.get().coordinates(), quad.get().weights());
  std::cout << "V on P3 quad with no filter:" << std::endl;
  std::cout << fe.Vandermonde() << std::endl;

  fe.fill_Vandermonde(quad.get().coordinates(), quad.get().weights(), true);
  std::cout << "V on P3 quad with filter:" << std::endl;
  std::cout << fe.Vandermonde() << std::endl;
}

// ----------------------------------------------------------------------------
