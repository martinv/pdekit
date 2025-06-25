/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE trilinos_crs_matrix_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

/// PDEKIT headers
#include "common/MPI/MPIEnv.hpp"
#include "linear_system/TpetraCrsMatrix.hpp"
#include "linear_system/TrilinosCrsMatrix.hpp"

using namespace pdekit;
using namespace pdekit::common;

// ----------------------------------------------------------------------------

struct TrilinosCrsMatrix_Fixture
{
  /// common setup for each test case
  TrilinosCrsMatrix_Fixture()
  {
    // uncomment if you want to use arguments to the test executable
    m_argc = boost::unit_test::framework::master_test_suite().argc;
    m_argv = boost::unit_test::framework::master_test_suite().argv;
  }

  /// common tear-down for each test case
  ~TrilinosCrsMatrix_Fixture()
  {
  }

  /// possibly common variables/functions used on the tests below
  int m_argc;
  char **m_argv;
};

BOOST_FIXTURE_TEST_SUITE(Trilinos_CRS_Matrix_TestSuite, TrilinosCrsMatrix_Fixture)

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(parallel_environment_init)
{
  common::mpi::MPIEnv::instance_type &mpi_env = common::mpi::MPIEnv::instance(m_argc, m_argv);
  BOOST_CHECK_EQUAL(mpi_env.is_initialized(), true);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_crs_matrix_utest)
{
#if PDEKIT_HAVE_TRILINOS

  // Create tridiagonal matrix with five rows
  const Uint nb_global_elem = 13;
  ls::TrilinosCrsMatrix crs_mat(nb_global_elem);

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {0.0, 0.0};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  indices = {3, 4};
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {0.0, 0.0, 0.0};
  indices = {0, 1, 2};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  indices = {1, 2, 3};
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock_structure();

  crs_mat.print_structure_to_file("sparsity.ps");

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_utest)
{
#if PDEKIT_HAVE_TRILINOS

  // Create tridiagonal matrix with five rows
  const Uint nb_global_elem = 5;
  ls::TpetraCrsMatrix<Real> crs_mat(nb_global_elem);

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {0.0, 0.0};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  indices = {3, 4};
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {0.0, 0.0, 0.0};
  indices = {0, 1, 2};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  indices = {1, 2, 3};
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock();

// crs_mat.print_structure_to_file("sparsity.ps");

// std::cout << crs_mat << std::endl;
#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_fixed_sparsity_from_graph_utest)
{
#if PDEKIT_HAVE_TRILINOS

  // Create tridiagonal matrix with five rows
  //  { 0, 1 }, { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4 }
  graph::Graph<Int> nonzero_positions(5);
  nonzero_positions.insert_edge(0, 0);
  nonzero_positions.insert_edge(0, 1);

  nonzero_positions.insert_edge(1, 0);
  nonzero_positions.insert_edge(1, 1);
  nonzero_positions.insert_edge(1, 2);

  nonzero_positions.insert_edge(2, 1);
  nonzero_positions.insert_edge(2, 2);
  nonzero_positions.insert_edge(2, 3);

  nonzero_positions.insert_edge(3, 2);
  nonzero_positions.insert_edge(3, 3);
  nonzero_positions.insert_edge(3, 4);

  nonzero_positions.insert_edge(4, 3);
  nonzero_positions.insert_edge(4, 4);

  /*
  ls::TpetraCrsMatrix<Real> crs_mat(nonzero_positions);
  */
  ls::TpetraCrsMatrix<Real> crs_mat;
  crs_mat.init(nonzero_positions);
  crs_mat.unlock();

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {0.0, 0.0};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  values  = {0.0, 0.1};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  values  = {4.3, 4.4};
  indices = {3, 4};
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {1.0, 1.1, 1.2};
  indices = {0, 1, 2};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  values  = {2.1, 2.2, 2.3};
  indices = {1, 2, 3};
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  values  = {3.2, 3.3, 3.4};
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock();
  // std::cout << crs_mat << std::endl;

  math::DenseConstVecView<Int> row_indices_view;
  math::DenseConstVecView<Real> row_values_view;
  crs_mat.row(0, row_indices_view, row_values_view);
  std::cout << "Indices:" << std::endl;
  std::cout << row_indices_view << std::endl;
  std::cout << "Values:" << std::endl;
  std::cout << row_values_view << std::endl;

  std::vector<Real> block_values(4);
  math::DenseMatView<Real> block_view(block_values.data(), 2, 2, 2);

  crs_mat.get_block(common::Range1D<Int>(2, 3), common::Range1D<Int>(2, 3), block_view);
  std::cout << block_view << std::endl;

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_fixed_sparsity_from_pattern_utest)
{
#if PDEKIT_HAVE_TRILINOS

  const std::vector<std::vector<Int>> sparsity = {{0, 1}, {0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 4}};
  math::MatrixSparsityPattern<Int> spattern(5, 5);
  spattern.build_sparsity(sparsity);

  ls::TpetraCrsMatrix<Real> crs_mat;
  crs_mat.init(spattern);
  crs_mat.unlock();

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {0.0, 0.0};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  values  = {0.0, 0.1};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  values  = {4.3, 4.4};
  indices = {3, 4};
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {1.0, 1.1, 1.2};
  indices = {0, 1, 2};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  values  = {2.1, 2.2, 2.3};
  indices = {1, 2, 3};
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  values  = {3.2, 3.3, 3.4};
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock();
  // std::cout << crs_mat << std::endl;

  math::DenseConstVecView<Int> row_indices_view;
  math::DenseConstVecView<Real> row_values_view;
  crs_mat.row(0, row_indices_view, row_values_view);
  std::cout << "Indices:" << std::endl;
  std::cout << row_indices_view << std::endl;
  std::cout << "Values:" << std::endl;
  std::cout << row_values_view << std::endl;

  std::vector<Real> block_values(4);
  math::DenseMatView<Real> block_view(block_values.data(), 2, 2, 2);

  crs_mat.get_block(common::Range1D<Int>(2, 3), common::Range1D<Int>(2, 3), block_view);
  std::cout << block_view << std::endl;

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_fixed_sparsity_from_block_pattern_utest)
{
#if PDEKIT_HAVE_TRILINOS

  std::unique_ptr<common::BlockArray<std::tuple<Int, Int>, Int>> nonzero_blocks(
      new common::BlockArray<std::tuple<Int, Int>, Int>());
  std::unique_ptr<std::vector<std::tuple<Int, Int>>> block_line_positions(
      new std::vector<std::tuple<Int, Int>>());

  const std::vector<std::tuple<Int, Int>> block_values = {{0, 3}, {4, 6}};
  const std::vector<std::tuple<Int, Int>> row_limits   = {{0, 3}, {4, 6}};

  nonzero_blocks->resize(4, 2);
  const common::ArrayView<const std::tuple<Int, Int>, _1D, Int> block_values_view0(
      block_values.data(), 1);
  nonzero_blocks->insert_block(0, block_values_view0);
  const common::ArrayView<const std::tuple<Int, Int>, _1D, Int> block_values_view1(
      block_values.data() + 1, 1);
  nonzero_blocks->insert_block(1, block_values_view1);

  // std::cout << (*nonzero_blocks) << std::endl;

#if 0
  const std::vector<std::vector<Int>> sparsity = {
    { 0, 1 }, { 0, 1, 2 }, { 1, 2, 3 }, { 2, 3, 4 }, { 3, 4 }
  };
  math::MatrixSparsityPattern<Int> spattern(5, 5);
  spattern.build_sparsity(sparsity);

  ls::TpetraCrsMatrix<Real> crs_mat;
  crs_mat.init(spattern);
  crs_mat.unlock();

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values = { 0.0, 0.0 };
  indices = { 0, 1 };
  crs_mat.insert_values_in_row(0, values, indices);

  values = { 0.0, 0.1 };
  indices = { 0, 1 };
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  values = { 4.3, 4.4 };
  indices = { 3, 4 };
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values = { 1.0, 1.1, 1.2 };
  indices = { 0, 1, 2 };
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  values = { 2.1, 2.2, 2.3 };
  indices = { 1, 2, 3 };
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = { 2, 3, 4 };
  values = { 3.2, 3.3, 3.4 };
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock();
  // std::cout << crs_mat << std::endl;

  math::DenseConstVecView<Int> row_indices_view;
  math::DenseConstVecView<Real> row_values_view;
  crs_mat.row(0, row_indices_view, row_values_view);
  std::cout << "Indices:" << std::endl;
  std::cout << row_indices_view << std::endl;
  std::cout << "Values:" << std::endl;
  std::cout << row_values_view << std::endl;

  std::vector<Real> block_values(4);
  math::DenseMatView<Real> block_view(block_values.data(), 2, 2, 2);

  crs_mat.get_block(std::make_tuple(2, 2), std::make_tuple(3, 3), block_view);
  std::cout << block_view << std::endl;
#endif

#endif
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(trilinos_tpetra_matrix_as_operator_utest)
{
#if PDEKIT_HAVE_TRILINOS

  const std::vector<std::vector<Int>> sparsity = {{0, 1}, {0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 4}};
  math::MatrixSparsityPattern<Int> spattern(5, 5);
  spattern.build_sparsity(sparsity);

  ls::TpetraCrsMatrix<Real> crs_mat;
  crs_mat.init(spattern);
  crs_mat.unlock();

  std::vector<Real> values;
  std::vector<Int> indices;

  // Insert first row
  values  = {0.0, 0.0};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  values  = {0.0, 0.1};
  indices = {0, 1};
  crs_mat.insert_values_in_row(0, values, indices);

  // Insert last row
  values  = {4.3, 4.4};
  indices = {3, 4};
  crs_mat.insert_values_in_row(4, values, indices);

  // Insert second row
  values  = {1.0, 1.1, 1.2};
  indices = {0, 1, 2};
  crs_mat.insert_values_in_row(1, values, indices);

  // Insert third row
  values  = {2.1, 2.2, 2.3};
  indices = {1, 2, 3};
  crs_mat.insert_values_in_row(2, values, indices);

  // Insert fourth row
  indices = {2, 3, 4};
  values  = {3.2, 3.3, 3.4};
  crs_mat.insert_values_in_row(3, values, indices);

  crs_mat.lock();
  // std::cout << crs_mat << std::endl;

  ls::TpetraMultiVector<Real> vec_X(crs_mat.domain_map());
  ls::TpetraMultiVector<Real> vec_Y(crs_mat.range_map());
  vec_X.insert_value(0, 0.0);
  vec_X.insert_value(1, 1.0);
  vec_X.insert_value(2, 2.0);
  vec_X.insert_value(3, 3.0);
  vec_X.insert_value(4, 4.0);

  vec_Y.insert_value(0, 0.0);
  vec_Y.insert_value(1, 0.0);
  vec_Y.insert_value(2, 0.0);
  vec_Y.insert_value(3, 0.0);
  vec_Y.insert_value(4, 0.0);

  // Compute y = 1.0 * y + 1.0 * A * x
  crs_mat.apply(vec_X, vec_Y, false, 1.0, 1.0);

  BOOST_CHECK_CLOSE(vec_Y.value(0), 0.1, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(1), 3.5, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(2), 13.4, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(3), 29.9, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(4), 30.5, 1.e-13);

  // Compute y = 0.0 * y + 1.0 * transpose(A) * x
  crs_mat.apply(vec_X, vec_Y, true);

  BOOST_CHECK_CLOSE(vec_Y.value(0), 1.0, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(1), 5.3, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(2), 15.2, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(3), 31.7, 1.e-13);
  BOOST_CHECK_CLOSE(vec_Y.value(4), 27.8, 1.e-13);

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
