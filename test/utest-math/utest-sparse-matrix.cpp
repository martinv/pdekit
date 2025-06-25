/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sparse_matrix_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/MatrixSparsityPattern.hpp"
#include "math/SparseDMat.hpp"

using namespace pdekit;
using namespace pdekit::math;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(crs_sparse_matrix_pattern_indexing_utest)
{
  std::vector<Uint> raw_row_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  std::unique_ptr<common::BlockArray<Uint, Uint>> compressed_rows(
      new common::BlockArray<Uint, Uint>());
  BOOST_CHECK_EQUAL(compressed_rows->nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(compressed_rows->size(), 0u);

  compressed_rows->reserve(9, 3);

  compressed_rows->create_back_block(3);
  compressed_rows->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_row_values.data(), 3));

  compressed_rows->create_back_block(3);
  compressed_rows->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_row_values.data() + 3, 3));

  compressed_rows->create_back_block(3);
  compressed_rows->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_row_values.data() + 6, 3));

  math::MatrixSparsityPattern<Uint, RowMajor> spattern;
  spattern.build_sparsity(std::move(compressed_rows));
  // spattern.print_sparsity();

  /*
  spattern.print_svg("sparsity.svg");
  spattern.print_vtu("sparsity.vtu");
  */

  // Check the first row
  // The result res0, res1 etc. is the global index of value searched for
  // in compressed_rows. So if we search for value 2 in the first row, the
  // result should be 1, because that's the position of 2 in raw_row_values
  const Uint res0 = spattern.idx_position(0, 1);
  const Uint res1 = spattern.idx_position(0, 2);
  const Uint res2 = spattern.idx_position(0, 3);

  BOOST_CHECK_EQUAL(res0, 0u);
  BOOST_CHECK_EQUAL(res1, 1u);
  BOOST_CHECK_EQUAL(res2, 2u);

  // Check the second row
  const Uint res3 = spattern.idx_position(1, 4);
  const Uint res4 = spattern.idx_position(1, 5);
  const Uint res5 = spattern.idx_position(1, 6);

  BOOST_CHECK_EQUAL(res3, 3u);
  BOOST_CHECK_EQUAL(res4, 4u);
  BOOST_CHECK_EQUAL(res5, 5u);

  // Check the third row
  const Uint res6 = spattern.idx_position(2, 7);
  const Uint res7 = spattern.idx_position(2, 8);
  const Uint res8 = spattern.idx_position(2, 9);

  BOOST_CHECK_EQUAL(res6, 6u);
  BOOST_CHECK_EQUAL(res7, 7u);
  BOOST_CHECK_EQUAL(res8, 8u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ccs_sparse_matrix_indexing_utest)
{
  std::vector<Uint> raw_col_values = {10, 20, 30, 40, 50, 60, 70, 80, 90};

  std::unique_ptr<common::BlockArray<Uint, Uint>> compressed_cols(
      new common::BlockArray<Uint, Uint>());
  BOOST_CHECK_EQUAL(compressed_cols->nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(compressed_cols->size(), 0u);

  compressed_cols->reserve(9, 3);

  compressed_cols->create_back_block(4);
  compressed_cols->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_col_values.data(), 4));

  compressed_cols->create_back_block(2);
  compressed_cols->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_col_values.data() + 4, 2));

  compressed_cols->create_back_block(3);
  compressed_cols->fill_last_block(
      common::ArrayView<const Uint, _1D, Uint>(raw_col_values.data() + 6, 3));

  math::MatrixSparsityPattern<Uint, ColumnMajor> spattern;
  spattern.build_sparsity(std::move(compressed_cols));
  // spattern.print_sparsity();

  // Check first column
  // The result res0, res1 etc. is the global index of value searched for
  // in compressed_cols. So if we search for value 20 in the first column, the
  // result should be 1, because that's the position of 20 in raw_col_values
  const Uint res0 = spattern.idx_position(10, 0);
  const Uint res1 = spattern.idx_position(20, 0);
  const Uint res2 = spattern.idx_position(30, 0);
  const Uint res3 = spattern.idx_position(40, 0);

  BOOST_CHECK_EQUAL(res0, 0u);
  BOOST_CHECK_EQUAL(res1, 1u);
  BOOST_CHECK_EQUAL(res2, 2u);
  BOOST_CHECK_EQUAL(res3, 3u);

  // Check second column
  const Uint res4 = spattern.idx_position(50, 1);
  const Uint res5 = spattern.idx_position(60, 1);

  BOOST_CHECK_EQUAL(res4, 4u);
  BOOST_CHECK_EQUAL(res5, 5u);

  // Check third column
  const Uint res6 = spattern.idx_position(70, 2);
  const Uint res7 = spattern.idx_position(80, 2);
  const Uint res8 = spattern.idx_position(90, 2);

  BOOST_CHECK_EQUAL(res6, 6u);
  BOOST_CHECK_EQUAL(res7, 7u);
  BOOST_CHECK_EQUAL(res8, 8u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(sparse_matrix_construct_utest)
{
  math::SparseDMat<Real> sparse_mat;

  std::vector<std::vector<Uint>> nonzeros;
  nonzeros.push_back({0, 1});
  nonzeros.push_back({0, 1, 2});
  nonzeros.push_back({1, 2, 3});
  nonzeros.push_back({2, 3});

  sparse_mat.resize(4, 4);
  sparse_mat.build_sparsity(nonzeros);

  Uint nnz = 0;
  for (const std::vector<Uint> &line : nonzeros)
  {
    nnz += line.size();
  }

  std::vector<Real> sparse_entries(sparse_mat.nb_nz());
  for (Uint i = 0; i < sparse_entries.size(); ++i)
  {
    sparse_entries[i] = 1.4 + i * 0.2;
  }

  const math::DenseConstVecView<Real> entries_view(sparse_entries.data(), sparse_entries.size());
  sparse_mat.fill_data(entries_view);

  sparse_mat.print();

  const std::tuple<const common::ArrayView<const Uint, _1D, Uint>,
                   const math::DenseConstVecView<Real>>
      const_line_data = sparse_mat.const_line_data(3);
  std::cout << "Indices on line 1 = " << std::get<0>(const_line_data) << std::endl;
  std::cout << "Values on line 1 = " << std::get<1>(const_line_data) << std::endl;

  std::tuple<const common::ArrayView<const Uint, _1D, Uint>, math::DenseVecView<Real>>
      non_const_line_data = sparse_mat.line_data(3);

  math::DenseVecView<Real> line = std::get<1>(non_const_line_data);
  line[0]                       = 100.0;

  sparse_mat.print();

  // auto raw_data = sparse_mat.raw_data();
}

// ----------------------------------------------------------------------------
