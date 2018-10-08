/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE singleton_test
#include <boost/test/unit_test.hpp>

#include <iostream>

#include "common/PDEKit.hpp"
#include "common/Singleton.hpp"

using namespace pdekit;
using namespace pdekit::common;

class Test
{
  public:
  void print()
  {
    std::cout << "This is test class" << std::endl;
  }
};

typedef Singleton<Test> SingletonTest;

BOOST_AUTO_TEST_CASE(singleton_utest)
{

  SingletonTest::instance_type &t = SingletonTest::instance();

  t.print();
}
