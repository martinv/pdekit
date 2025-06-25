/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE environment_utest
#include <boost/bind.hpp>
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/Core.hpp"
#include "common/MPI/MPIEnv.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

struct Core_Fixture
{
  /// common setup for each test case
  Core_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~Core_Fixture()
  {
  }

  /// possibly common functions used on the tests below
  int m_argc;
  char **m_argv;
};

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(Core_TestSuite, Core_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(init_environment)
{
  std::cout << m_argc << std::endl;
  std::cout << m_argv[0] << std::endl;

  common::Core::instance_type &environment = common::Core::instance(m_argc, m_argv);
  // common::MPI::MPIEnv::instance_type& par_env =
  // common::MPI::MPIEnv::instance(m_argc,m_argv);

  // par_env.finalize();

  environment.set_nb_threads(2);
  BOOST_CHECK_EQUAL(environment.nb_threads(), 2u);

  environment.finalize();
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()

// ----------------------------------------------------------------------------
