/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE trilinos_ls_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "linear_system/LSTrilinos.hpp"

#include "linear_system/LSTpetra.hpp"

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

BOOST_AUTO_TEST_CASE(trilinos_ls_utest_parenv_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_ls_utest)
{
#if PDEKIT_HAVE_TRILINOS

  const Uint nb_rows = 5;

  std::shared_ptr<TrilinosCrsMatrix> matrix(new TrilinosCrsMatrix(nb_rows));

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {2.0, 1.0};
  indices = {0, 1};
  matrix->insert_values_in_row(0, values, indices);

  // Insert last row
  values  = {1.0, 2.0};
  indices = {3, 4};
  matrix->insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {1.0, 2.0, 1.0};
  indices = {0, 1, 2};
  matrix->insert_values_in_row(1, values, indices);

  // Insert third row
  indices = {1, 2, 3};
  matrix->insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  matrix->insert_values_in_row(3, values, indices);

  matrix->lock_structure();
  // std::cout << *matrix << std::endl;

  std::shared_ptr<TrilinosMultiVector> b(new TrilinosMultiVector(*matrix->map()));
  BOOST_CHECK_EQUAL(b->size(), nb_rows);

  (*b)(0) = 0.4;
  (*b)(1) = 0.8;
  (*b)(2) = 1.2;
  (*b)(3) = 1.6;
  (*b)(4) = 1.4;

  // std::cout << *b << std::endl;

  std::shared_ptr<TrilinosMultiVector> x(new TrilinosMultiVector(*matrix->map()));
  BOOST_CHECK_EQUAL(x->size(), nb_rows);

  LSTrilinos lin_system;
  lin_system.configure(matrix, b, x);

  x->fill(0.0);
  lin_system.solve();

  const Real tol = 1.e-12;

  BOOST_CHECK_CLOSE((*x)(0), 0.1, tol);
  BOOST_CHECK_CLOSE((*x)(1), 0.2, tol);
  BOOST_CHECK_CLOSE((*x)(2), 0.3, tol);
  BOOST_CHECK_CLOSE((*x)(3), 0.4, tol);
  BOOST_CHECK_CLOSE((*x)(4), 0.5, tol);

  // std::cout << *x << std::endl;

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(tpetra_ls_utest)
{
#if PDEKIT_HAVE_TRILINOS

  const Uint nb_rows = 5;

  std::shared_ptr<TpetraCrsMatrix<Real>> matrix(new TpetraCrsMatrix<Real>(nb_rows));

  std::vector<Real> values;
  std::vector<Int> indices;

  // Fill system matrix
  // Insert first row
  values  = {2.0, 1.0};
  indices = {0, 1};
  matrix->insert_values_in_row(0, values, indices);

  // Insert last row
  values  = {1.0, 2.0};
  indices = {3, 4};
  matrix->insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {1.0, 2.0, 1.0};
  indices = {0, 1, 2};
  matrix->insert_values_in_row(1, values, indices);

  // Insert third row
  indices = {1, 2, 3};
  matrix->insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  matrix->insert_values_in_row(3, values, indices);

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

  LSTpetra<Real> lin_system;
  // lin_system.configure(matrix, b, x, true, false);
  lin_system.initialize_solver(matrix, b, x, false);
  lin_system.update_after_mat_values_change(matrix, b, x, false);

  x->fill(0.0);
  lin_system.solve();

  const Real tol = 1.e-12;

  vec_column = 0;

  BOOST_CHECK_CLOSE((*x).value(0, vec_column), 0.1, tol);
  BOOST_CHECK_CLOSE((*x).value(1, vec_column), 0.2, tol);
  BOOST_CHECK_CLOSE((*x).value(2, vec_column), 0.3, tol);
  BOOST_CHECK_CLOSE((*x).value(3, vec_column), 0.4, tol);
  BOOST_CHECK_CLOSE((*x).value(4, vec_column), 0.5, tol);

  vec_column = 1;

  BOOST_CHECK_CLOSE((*x).value(0, vec_column), 0.2, tol);
  BOOST_CHECK_CLOSE((*x).value(1, vec_column), 0.4, tol);
  BOOST_CHECK_CLOSE((*x).value(2, vec_column), 0.6, tol);
  BOOST_CHECK_CLOSE((*x).value(3, vec_column), 0.8, tol);
  BOOST_CHECK_CLOSE((*x).value(4, vec_column), 1.0, tol);

  //  std::cout << *x << std::endl;

  x->fill(1.0);
  std::shared_ptr<TpetraMultiVector<Real>> residual(
      new TpetraMultiVector<Real>(matrix->map(), nb_rhs_columns));
  BOOST_CHECK_EQUAL(residual->size(), nb_rows);

  lin_system.residual(*residual);

  vec_column = 0;
  BOOST_CHECK_CLOSE((*residual).value(0, vec_column), -2.6, tol);
  BOOST_CHECK_CLOSE((*residual).value(1, vec_column), -3.2, tol);
  BOOST_CHECK_CLOSE((*residual).value(2, vec_column), -2.8, tol);
  BOOST_CHECK_CLOSE((*residual).value(3, vec_column), -2.4, tol);
  BOOST_CHECK_CLOSE((*residual).value(4, vec_column), -1.6, tol);

  vec_column = 1;
  BOOST_CHECK_CLOSE((*residual).value(0, vec_column), -2.2, tol);
  BOOST_CHECK_CLOSE((*residual).value(1, vec_column), -2.4, tol);
  BOOST_CHECK_CLOSE((*residual).value(2, vec_column), -1.6, tol);
  BOOST_CHECK_CLOSE((*residual).value(3, vec_column), -0.8, tol);
  BOOST_CHECK_CLOSE((*residual).value(4, vec_column), -0.2, tol);

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_ls_utest_parenv_finalize)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance();
  mpi_env.finalize();
  BOOST_CHECK_EQUAL(mpi_env.is_finalized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END()
