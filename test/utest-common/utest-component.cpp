/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE "Test module for Component"
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/Component.hpp"
#include "common/StringUtils.hpp"

using namespace pdekit;
using namespace pdekit::common;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

template <Uint N>
class TestComponent : public Component
{
  public:
  /// Constructor
  TestComponent(const std::string &name);

  /// Destructor
  ~TestComponent() override;

  /// Return type name of this component
  std::string derived_type_name() const override;

  private:
};

// ----------------------------------------------------------------------------

template <Uint N>
TestComponent<N>::TestComponent(const std::string &name) : Component(name)
{
}

template <Uint N>
TestComponent<N>::~TestComponent()
{
}

template <Uint N>
std::string TestComponent<N>::derived_type_name() const
{
  return "TestComponent<" + StringUtils::to_string(N) + ">";
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Component_TestSuite)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(component_constructors)
{
  std::shared_ptr<Component> root = make_component<TestComponent<0>>("root");

  BOOST_CHECK_EQUAL(root->name(), "root");
  BOOST_CHECK_EQUAL(root->uri().base_path().string(), "cpath:/");
  BOOST_CHECK_EQUAL(root->uri().string(), "cpath:/");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(adding_child_components)
{
  std::shared_ptr<Component> root = make_component<TestComponent<0>>("root");
  std::shared_ptr<Component> dir1 = make_component<TestComponent<1>>("dir1");
  std::shared_ptr<Component> dir2 = make_component<TestComponent<2>>("dir2");

  root->add_component(dir1);
  dir1->add_component(dir2);

  BOOST_CHECK_EQUAL(root->uri().string(), "cpath:/");
  BOOST_CHECK_EQUAL(dir1->uri().string(), "cpath:/dir1");
  BOOST_CHECK_EQUAL(dir2->uri().string(), "cpath:/dir1/dir2");

  // std::cout << root->tree() << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(access_component)
{
  std::shared_ptr<Component> root = make_component<TestComponent<0>>("root");

  std::shared_ptr<Component> dir1  = make_component<TestComponent<1>>("dir1");
  std::shared_ptr<Component> dir2  = make_component<TestComponent<2>>("dir2");
  std::shared_ptr<Component> dir21 = make_component<TestComponent<3>>("dir21");
  std::shared_ptr<Component> dir22 = make_component<TestComponent<3>>("dir22");

  // add child components to root
  root->add_component(dir1);
  dir1->add_component(dir2);
  dir2->add_component(dir21);
  dir2->add_component(dir22);

  // test relative & complete path
  URI p0("cpath:../dir21");
  std::shared_ptr<Component> cp0 = dir22->access_component(p0);
  BOOST_CHECK_EQUAL(cp0->uri().string(), "cpath:/dir1/dir2/dir21");

  // test relative & complete path
  URI p1("cpath:/dir1");
  std::shared_ptr<Component> cp1 = dir22->access_component(p1);
  BOOST_CHECK_EQUAL(cp1->uri().string(), "cpath:/dir1");
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
