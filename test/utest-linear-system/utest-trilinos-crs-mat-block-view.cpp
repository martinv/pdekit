/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tpetra_crs_matrix_block_view_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "linear_system/TpetraCrsMatrixDiagBlocks.hpp"

using namespace pdekit;
using namespace pdekit::common;

// ----------------------------------------------------------------------------

struct TpetraCrsMatrix_Fixture
{
  /// common setup for each test case
  TpetraCrsMatrix_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~TpetraCrsMatrix_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

BOOST_FIXTURE_TEST_SUITE(Tpetra_CRS_Matrix_TestSuite, TpetraCrsMatrix_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(parallel_environment_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_utest)
{
#if PDEKIT_HAVE_TRILINOS

  // Create tridiagonal matrix with five rows
  const Uint nb_global_elem = 16;
  ls::TpetraCrsMatrix<Real> crs_mat(nb_global_elem);

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert rows 1 and 2 - block 1
  indices = {0, 1};
  values  = {1.0, 2.0};
  crs_mat.insert_values_in_row(0, values, indices);
  values = {3.0, 4.0};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert rows 3 and 4 - block 2
  indices = {2, 3};
  values  = {5.0, 6.0};
  crs_mat.insert_values_in_row(2, values, indices);
  values = {7.0, 8.0};
  crs_mat.insert_values_in_row(3, values, indices);

  // Insert rows 5 and 6 - block 3
  indices = {4, 5};
  values  = {9.0, 10.0};
  crs_mat.insert_values_in_row(4, values, indices);
  values = {11.0, 12.0};
  crs_mat.insert_values_in_row(5, values, indices);

  indices = {6, 7};
  values  = {13.0, 14.0};
  crs_mat.insert_values_in_row(6, values, indices);
  values = {15.0, 16.0};
  crs_mat.insert_values_in_row(7, values, indices);

  crs_mat.lock();

  // crs_mat.print_structure_to_file("sparsity.ps");

  // std::cout << crs_mat << std::endl;

  ls::TpetraCrsMatrixDiagBlocks<Real> blocks1, blocks2;

  blocks1.set_block_sizes({2, 2, 2, 2});
  std::vector<Real> block_values;

  blocks1.get_block_data(crs_mat, 0, values);
  BOOST_CHECK_EQUAL(values[0], 1.0);
  BOOST_CHECK_EQUAL(values[1], 2.0);
  BOOST_CHECK_EQUAL(values[2], 3.0);
  BOOST_CHECK_EQUAL(values[3], 4.0);

  blocks1.get_block_data(crs_mat, 1, values);
  BOOST_CHECK_EQUAL(values[0], 5.0);
  BOOST_CHECK_EQUAL(values[1], 6.0);
  BOOST_CHECK_EQUAL(values[2], 7.0);
  BOOST_CHECK_EQUAL(values[3], 8.0);

  blocks1.get_block_data(crs_mat, 2, values);
  BOOST_CHECK_EQUAL(values[0], 9.0);
  BOOST_CHECK_EQUAL(values[1], 10.0);
  BOOST_CHECK_EQUAL(values[2], 11.0);
  BOOST_CHECK_EQUAL(values[3], 12.0);

  blocks1.get_block_data(crs_mat, 3, values);
  BOOST_CHECK_EQUAL(values[0], 13.0);
  BOOST_CHECK_EQUAL(values[1], 14.0);
  BOOST_CHECK_EQUAL(values[2], 15.0);
  BOOST_CHECK_EQUAL(values[3], 16.0);

  blocks1 = blocks2;
#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(parallel_environment_finalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
