/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE mpi_send_wrapper_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "common/MPI/MPIRecv.hpp"
#include "common/MPI/MPISend.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

struct MPI_Send_Fixture
{
  /// common setup for each test case
  MPI_Send_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~MPI_Send_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

// ----------------------------------------------------------------------------

template <typename T>
void check_vector_mpi_send_recv(const std::vector<T> &v_in)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();

  BOOST_CHECK_EQUAL(mpi_env.comm_size(), 2u);

  const auto my_rank      = mpi_env.rank();
  const auto communicator = mpi_env.comm();

  if (my_rank == 0)
  {
    common::mpi::send(v_in, 1, 10, communicator);
  }
  else
  {
    std::vector<T> v_out;
    common::mpi::recv(v_out, 0, 10, communicator);
    std::cout << "v_out = [ ";
    for (const auto val : v_out)
    {
      std::cout << val << " ";
    }
    std::cout << "]" << std::endl;

    BOOST_CHECK_EQUAL(v_out.size(), 5u);

    for (Uint i = 0; i < v_in.size(); ++i)
    {
      BOOST_CHECK_EQUAL(v_in[i], v_out[i]);
    }
  }
}

// ----------------------------------------------------------------------------

BOOST_FIXTURE_TEST_SUITE(MPI_Send_TestSuite, MPI_Send_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mpi_send_utest_env_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mpi_int_vector_send_recv_utest)
{
  const std::vector<Int> v_in = {101, 201, 301, 401, 501};
  check_vector_mpi_send_recv(v_in);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mpi_real_vector_send_recv_utest)
{
  const std::vector<Real> v_in = {1.1, 2.1, 3.1, 4.1, 5.1};
  check_vector_mpi_send_recv(v_in);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(mpi_send_utest_env_finalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
