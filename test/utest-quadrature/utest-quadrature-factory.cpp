/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE quadrature_factory_test
#include <boost/test/unit_test.hpp>

#include <cmath>
#include <iostream>

#include "common/PDEKit.hpp"
#include "mesh/point_set/StdPointSet.hpp"
#include "mesh/point_set/StdPointSetFactory.hpp"

#include <boost/shared_ptr.hpp>

using namespace pdekit;
using namespace pdekit::common;

BOOST_AUTO_TEST_CASE(quadrature_factory_utest)
{

  // FactoryT<QuadratureBase> fact;

  std::cout << "Asking for an instance:" << std::endl;
  StdPointSetFactory::instance_type &factory = StdPointSetFactory::instance();
  std::cout << " ... got the first instance" << std::endl;
  // FactoryT<ShapeFunction>& factory = FactoryT<ShapeFunction>::instance();
  // LagrangeP1Triag2D sf1;
  factory.print_builders();

  std::cout << "Asking for another instance:" << std::endl;
  StdPointSetFactory::instance_type &factory2 = StdPointSetFactory::instance();
  std::cout << "... got the second instance" << std::endl;
  // FactoryT<ShapeFunction>& factory2 = FactoryT<ShapeFunction>::instance();
  // LagrangeP2Triag2D sf2;
  factory2.print_builders();

  // FactoryT<QuadratureBase>::BuilderBase* builder =
  // factory.builder("Triag-P1");
  // QuadBuilderBase* builder = factory.builder("Triag-P1-Gauss");
  // SFFactory::BuilderBase* builder2 = factory.builder("LagrangeP2Triag2D");

  const StdPointSetFactory::instance_type::const_product_base_ptr quad =
      factory.create(mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Gauss));

  // builder->release_memory(quad);

  FactoryPool::instance_type &f_p = FactoryPool::instance();

  f_p.list_factories();

  const StdPointSetFactory::instance_type::const_product_base_ptr quad2 =
      factory.create(mesh::PointSetTag(ElemShape::Triag, P1, PointSetID::Gauss));
  std::cout << quad2->name() << std::endl;

  // delete quad2;

  const StdPointSetFactory::instance_type::const_product_base_ptr quad3 =
      factory.create(mesh::PointSetTag(ElemShape::Quad, P3, PointSetID::Gauss));

  std::cout << "Sum of weights for quadrature " << quad3->name() << " = ";

  math::DenseDVec<Real> w;

  quad3->weights(w);

  Real sum = 0.0;

  for (Uint q = 0; q < quad3->size(); ++q)
  {
    sum += w[q];
  }

  std::cout << sum << std::endl;
}
