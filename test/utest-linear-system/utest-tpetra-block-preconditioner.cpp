/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tpetra_block_preconditioner_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "linear_system/LSTpetra.hpp"
#include "linear_system/LSTrilinos.hpp"
#include "linear_system/preconditioner/TpetraBlockDiagPC.hpp"

using namespace pdekit;

#if PDEKIT_HAVE_TRILINOS
using namespace pdekit::ls;
#endif

// ----------------------------------------------------------------------------

struct Trilinos_LS_Fixture
{
  /// common setup for each test case
  Trilinos_LS_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~Trilinos_LS_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

BOOST_FIXTURE_TEST_SUITE(Trilinos_LS_TestSuite, Trilinos_LS_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tpetra_block_preconditioner_utest_parenv_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tpetra_block_preconditioner_utest)
{
#if PDEKIT_HAVE_TRILINOS

  const Uint nb_rows = 5;

  std::shared_ptr<TpetraCrsMatrix<Real>> matrix(new TpetraCrsMatrix<Real>(nb_rows));

  std::vector<Real> values;
  std::vector<Int> indices;

  // Fill system matrix
  //
  // | 1 2  0  0  0 |
  // | 3 4  0  0  0 |
  // | 0 0  5  6  7 |
  // | 0 0  8  9 10 |
  // | 0 0 11 12  1 |

  // Insert first row
  values  = {1.0, 2.0};
  indices = {0, 1};
  matrix->insert_values_in_row(0, values, indices);

  // Insert second row
  values = {3.0, 4.0};
  matrix->insert_values_in_row(1, values, indices);

  // Insert third row
  values  = {5.0, 6.0, 7.0};
  indices = {2, 3, 4};
  matrix->insert_values_in_row(2, values, indices);

  // Insert fourth row
  values = {8.0, 9.0, 10.0};
  matrix->insert_values_in_row(3, values, indices);

  // Insert fifth row
  values = {11.0, 12.0, 1.0};
  matrix->insert_values_in_row(4, values, indices);

  matrix->lock();
  // std::cout << *matrix << std::endl;

  // The rhs vector b has 2 columns - we will solve the
  // the system twice for different right-hand side
  const Uint nb_rhs_columns = 2;

  std::shared_ptr<TpetraMultiVector<Real>> b(
      new TpetraMultiVector<Real>(matrix->map(), nb_rhs_columns));
  BOOST_CHECK_EQUAL(b->size(), nb_rows);

  // Fill the first right-hand side
  Uint vec_column = 0;

  (*b).insert_value(0, 0.4, vec_column);
  (*b).insert_value(1, 0.8, vec_column);
  (*b).insert_value(2, 1.2, vec_column);
  (*b).insert_value(3, 1.6, vec_column);
  (*b).insert_value(4, 1.4, vec_column);

  // Fill the second right-hand side
  vec_column = 1;

  (*b).insert_value(0, 0.8, vec_column);
  (*b).insert_value(1, 1.6, vec_column);
  (*b).insert_value(2, 2.4, vec_column);
  (*b).insert_value(3, 3.2, vec_column);
  (*b).insert_value(4, 2.8, vec_column);

  // std::cout << *b << std::endl;

  // The solution vector also has 2 columns corresponding
  // to the two different right-hand sides
  std::shared_ptr<TpetraMultiVector<Real>> x(
      new TpetraMultiVector<Real>(matrix->map(), nb_rhs_columns));
  BOOST_CHECK_EQUAL(x->size(), nb_rows);

  x->fill(0.0);

  //  std::cout << *x << std::endl;

#endif

  std::unique_ptr<std::vector<common::Range1D<Int>>> blocks(
      new std::vector<common::Range1D<Int>>());
  blocks->resize(2);
  (*blocks)[0] = common::Range1D<Int>(0, 1);
  (*blocks)[1] = common::Range1D<Int>(2, 4);

  std::shared_ptr<TpetraBlockDiagPC<Real>> preconditioner(new TpetraBlockDiagPC<Real>());

  preconditioner->create("", matrix, std::move(blocks));
  preconditioner->compute(matrix);
  // preconditioner->print();
  preconditioner->apply(*b, *x);

  const Real tol = 1.e-10;

  vec_column = 0;

  BOOST_CHECK_CLOSE((*x).value(0, vec_column), 0.0, tol);
  BOOST_CHECK_CLOSE((*x).value(1, vec_column), 0.2, tol);
  BOOST_CHECK_CLOSE((*x).value(2, vec_column), -0.35, tol);
  BOOST_CHECK_CLOSE((*x).value(3, vec_column), 0.433333333333, tol);
  BOOST_CHECK_CLOSE((*x).value(4, vec_column), 0.05, tol);

  vec_column = 1;

  BOOST_CHECK_CLOSE((*x).value(0, vec_column), 0.0, tol);
  BOOST_CHECK_CLOSE((*x).value(1, vec_column), 0.4, tol);
  BOOST_CHECK_CLOSE((*x).value(2, vec_column), -0.7, tol);
  BOOST_CHECK_CLOSE((*x).value(3, vec_column), 0.8666666666666, tol);
  BOOST_CHECK_CLOSE((*x).value(4, vec_column), 0.1, tol);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tpetra_block_preconditioner_utest_parenv_finalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
