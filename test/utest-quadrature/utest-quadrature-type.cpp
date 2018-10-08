/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_type_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/point_set/StdPointSet.hpp"

using namespace pdekit;

BOOST_AUTO_TEST_CASE(quadrature_type_utest)
{
  mesh::PointSetTag quad(ElemShape::Quad, P6, PointSetID::Gauss);

  std::string quad_desc;

  quad_desc = quad.as_string();
  std::cout << quad_desc << std::endl;

  quad_desc = "Tetra-P7-Lobatto";

  quad = mesh::PointSetTag::string_to_tag(quad_desc);

  std::cout << "Quadrature type:" << std::endl;
  std::cout << "[" << quad.as_string() << "]" << std::endl;

  // --------------------------------------------------------------------------

  std::cout << "The id of Lobatto quadrature for P1 Triag  = ";
  std::cout << mesh::PointSetTag::fields_to_string(ElemShape::Triag, P1, PointSetID::GaussLobatto)
            << std::endl;

  mesh::PointSetTag quad2(ElemShape::Triag, P1, PointSetID::GaussLobatto);

  ElemShape key1;
  Uint key2;
  PointSetID key3;
  quad2.decompose_into_fields(quad2, key1, key2, key3);

  // --------------------------------------------------------------------------
}
