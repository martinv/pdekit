/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE factory_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/FactoryT.hpp"
#include "common/PDEKit.hpp"
#include "common/Singleton.hpp"

#include <boost/shared_ptr.hpp>

using namespace pdekit;

// ----------------------------------------------------------------------------

struct FactoryFixture
{

  struct MyEntityBase
  {
    static const std::string type_name()
    {
      return "MyEntityBase";
    }

    virtual void print() const = 0;
  };

  struct MyEntity1 : public MyEntityBase
  {
    MyEntity1()
    {
    }
    void print() const override
    {
      std::cout << "MyEntity1" << std::endl;
    }
  };

  struct MyEntity2 : public MyEntityBase
  {
    MyEntity2()
    {
    }
    void print() const override
    {
      std::cout << "MyEntity2" << std::endl;
    }
  };

  struct MyEntity3 : public MyEntityBase
  {
    MyEntity3() : m_val1(0u), m_val2(0.0)
    {
    }

    MyEntity3(std::tuple<Uint, Real> const &params)
        : m_val1(std::get<0>(params)), m_val2(std::get<1>(params))
    {
    }

    void print() const override
    {
      std::cout << "MyEntity3" << std::endl;
      std::cout << "Value1 = " << m_val1 << std::endl;
      std::cout << "Value2 = " << m_val2 << std::endl;
    }

    Uint m_val1;
    Real m_val2;
  };

  struct MyFactorySetup
  {
    static void setup_instance(common::FactoryT<MyEntityBase, std::string> &factory);
  };

}; // End of fixture

void FactoryFixture::MyFactorySetup::setup_instance(
    common::FactoryT<MyEntityBase, std::string> &factory)
{
  factory.register_builder<MyEntity1>("MyEntity1");
  factory.register_builder<MyEntity2>("MyEntity2");
  factory.register_builder<MyEntity3>("MyEntity3a");
  factory.register_builder<MyEntity3, const Uint, const Real>("MyEntity3b", 11U, 8.3);

  Uint v1 = 12U;
  Real v2 = 9.3;
  factory.register_builder<MyEntity3, const Uint, const Real>("MyEntity3c", v1, v2);
}

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(Factory_TestSuite, FactoryFixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(factory_builders_registration_utest)
{

  typedef common::Singleton<common::FactoryT<MyEntityBase, std::string>, MyFactorySetup> MyFactory;

  MyFactory::instance_type &factory = MyFactory::instance();

  const std::vector<std::string> keys = {"MyEntity1", "MyEntity2", "MyEntity3a", "MyEntity3b",
                                         "MyEntity3c"};

  for (Uint k = 0; k < keys.size(); ++k)
  {
    const MyFactory::instance_type::const_product_base_ptr entity = factory.create(keys[k]);
    entity->print();
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
