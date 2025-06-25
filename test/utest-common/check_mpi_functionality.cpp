#define DONT_USE_BOOST_TEST 1
#if DONT_USE_BOOST_TEST
#include <iostream>

#include "common/MPI/MPIEnv.hpp"
#include "common/PDEKit.hpp"

using namespace pdekit;

int main(int argc, char *argv[])
{
  // MPI_Init(&argc,&argv);

  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(argc, argv);

  /*
  // Another option to initialize:
  common::mpi::MPIEnv::instance_type &mpi_env =
  common::mpi::MPIEnv::instance(); mpi_env.init(argc, argv);
  */

  const Uint rank = mpi_env.rank();

  if (rank == 0)
  {
    if (mpi_env.is_initialized())
    {
      std::cout << "MPI environment is properly initialized" << std::endl;
      std::cout << "MPI world communicator size = " << mpi_env.comm_size() << std::endl;
    }
    else
    {
      std::cerr << "Didn't manage to initialize the mpi envirnoment properly!" << std::endl;
    }
  }

  // MPI_Finalize();
  mpi_env.finalize();

  if (rank == 0)
  {
    if (mpi_env.is_finalized())
    {
      std::cout << "MPI environment was properly finalized" << std::endl;
    }
    else
    {
      std::cerr << "Didn't manage to finalize the mpi envirnoment properly!" << std::endl;
    }
  }

  return 0;
}

#else

/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mpi_functionality_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"

using namespace pdekit;

// ============================================================================

struct MPI_Functionality_Fixture
{
  /// common setup for each test case
  MPI_Functionality_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~MPI_Functionality_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

BOOST_FIXTURE_TEST_SUITE(MPI_Functionality_TestSuite, MPI_Functionality_Fixture)

// ============================================================================

BOOST_AUTO_TEST_CASE(mpi_functionality_parenv_init)
{
  // common::mpi::MPIEnv::instance_type& mpi_env =
  // common::mpi::MPIEnv::instance(m_argc,m_argv);
  // BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);

  MPI_Init(&m_argc, &m_argv);
}

// ============================================================================

BOOST_AUTO_TEST_CASE(mpi_functionality_parenv_finalize)
{
  // common::mpi::MPIEnv::instance_type& mpi_env =
  // common::mpi::MPIEnv::instance();
  // mpi_env.finalize();
  // BOOST_CHECK_EQUAL( mpi_env.is_finalized(), true );

  MPI_Finalize();
}

// ============================================================================

BOOST_AUTO_TEST_SUITE_END()
#endif

#if 0
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mpi_functionality
#include <boost/test/unit_test.hpp>

#include "mpi.h"

// ============================================================================

struct MPI_Functionality_Fixture
{
  MPI_Functionality_Fixture()
  {
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  int m_argc;
  char** m_argv;
};

BOOST_FIXTURE_TEST_SUITE( MPI_Functionality_TestSuite, MPI_Functionality_Fixture )

// ============================================================================

BOOST_AUTO_TEST_CASE( mpi_functionality_parenv_init )
{
  MPI_Init(&m_argc,&m_argv);
}

// ============================================================================

BOOST_AUTO_TEST_CASE( mpi_functionality_parenv_finalize )
{
  MPI_Finalize();
}

// ============================================================================

BOOST_AUTO_TEST_SUITE_END()
#endif
