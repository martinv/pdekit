/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE block_sparse_matrix_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

/// PDEKIT headers
#include "math/BlockMatrixSparsityPattern.hpp"

using namespace pdekit;
using namespace pdekit::math;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(crs_block_sparse_matrix_pattern_build_utest)
{
  std::vector<std::vector<common::Range1D<Int>>> raw_row_blocks(3);
  raw_row_blocks[0] = {{31, 70}, {0, 30}};
  raw_row_blocks[1] = {{31, 70}, {71, 110}};
  raw_row_blocks[2] = {{71, 110}};

  std::vector<std::vector<Int>> block_ids(3);
  block_ids[0] = {1, 0};
  block_ids[1] = {1, 2};
  block_ids[2] = {2};

  std::vector<common::Range1D<Int>> block_line_limits = {
      common::Range1D<Int>(0, 40), common::Range1D<Int>(41, 90), common::Range1D<Int>(91, 140)};

  math::BlockMatrixSparsityPattern<Int, RowMajor> spattern;
  spattern.build_sparsity(raw_row_blocks, block_ids, block_line_limits);

  const Uint nb_nz = 41 * (31 + 40) + 50 * (40 + 40) + 50 * (40);

  BOOST_CHECK_EQUAL(spattern.nb_nz(), nb_nz);
  BOOST_CHECK_EQUAL(spattern.nb_lines(), 3u);

  // ---------

  const Uint r0 = 24;
  const Uint c0 = 20;

  std::tuple<Uint, Uint> idx_pos0 = spattern.idx_position(r0, c0);
  BOOST_CHECK_EQUAL(std::get<0>(idx_pos0), 0u);
  BOOST_CHECK_EQUAL(std::get<1>(idx_pos0), 0u);

  // ---------

  const Uint r1 = 14;
  const Uint c1 = 70;

  std::tuple<Uint, Uint> idx_pos1 = spattern.idx_position(r1, c1);
  BOOST_CHECK_EQUAL(std::get<0>(idx_pos1), 0u);
  BOOST_CHECK_EQUAL(std::get<1>(idx_pos1), 1u);

  // ---------

  const Uint r2 = 45;
  const Uint c2 = 55;

  std::tuple<Uint, Uint> idx_pos2 = spattern.idx_position(r2, c2);
  BOOST_CHECK_EQUAL(std::get<0>(idx_pos2), 1u);
  BOOST_CHECK_EQUAL(std::get<1>(idx_pos2), 1u);

  // ---------

  const Uint r3 = 45;
  const Uint c3 = 99;

  std::tuple<Uint, Uint> idx_pos3 = spattern.idx_position(r3, c3);
  BOOST_CHECK_EQUAL(std::get<0>(idx_pos3), 1u);
  BOOST_CHECK_EQUAL(std::get<1>(idx_pos3), 2u);

  // ---------

  const Uint r4 = 140;
  const Uint c4 = 110;

  std::tuple<Uint, Uint> idx_pos4 = spattern.idx_position(r4, c4);
  BOOST_CHECK_EQUAL(std::get<0>(idx_pos4), 2u);
  BOOST_CHECK_EQUAL(std::get<1>(idx_pos4), 2u);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(crs_block_sparse_matrix_pattern_indexing_utest)
{
  std::vector<common::Range1D<Uint>> raw_row_blocks = {
      {0, 30}, {31, 70}, {31, 70}, {71, 110}, {71, 110}};

  std::unique_ptr<common::BlockMultiArray<common::Range1D<Uint>, Uint>> compressed_blk_rows(
      new common::BlockMultiArray<common::Range1D<Uint>, Uint>());
  BOOST_CHECK_EQUAL(compressed_blk_rows->nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(compressed_blk_rows->size(), 0u);

  compressed_blk_rows->reserve(5, 3);

  compressed_blk_rows->create_back_block(2);
  compressed_blk_rows->fill_last_block<0>(
      common::ArrayView<const common::Range1D<Uint>, _1D, Uint>(raw_row_blocks.data(), 2));

  compressed_blk_rows->create_back_block(2);
  compressed_blk_rows->fill_last_block<0>(
      common::ArrayView<const common::Range1D<Uint>, _1D, Uint>(raw_row_blocks.data() + 2, 2));

  compressed_blk_rows->create_back_block(1);
  compressed_blk_rows->fill_last_block<0>(
      common::ArrayView<const common::Range1D<Uint>, _1D, Uint>(raw_row_blocks.data() + 4, 1));

  std::vector<Uint> block_ids = {0, 1, 1, 2, 2};
  compressed_blk_rows->insert_block<1>(
      0, common::ArrayView<const Uint, _1D, Uint>(block_ids.data(), 2));
  compressed_blk_rows->insert_block<1>(
      1, common::ArrayView<const Uint, _1D, Uint>(block_ids.data() + 2, 2));
  compressed_blk_rows->insert_block<1>(
      2, common::ArrayView<const Uint, _1D, Uint>(block_ids.data() + 4, 1));

  std::unique_ptr<std::vector<common::Range1D<Uint>>> row_positions(
      new std::vector<common::Range1D<Uint>>());
  row_positions->push_back(common::Range1D<Uint>(0, 40));
  row_positions->push_back(common::Range1D<Uint>(41, 90));
  row_positions->push_back(common::Range1D<Uint>(91, 140));

  math::BlockMatrixSparsityPattern<Uint, RowMajor> spattern;
  spattern.build_sparsity(std::move(compressed_blk_rows), std::move(row_positions));

  // Check the first row
  // The result res0, res1 etc. is the global index of value searched for
  // in compressed_rows. So if we search for value 2 in the first row, the
  // result should be 1, because that's the position of 2 in raw_row_values
  const std::tuple<Uint, Uint> res0 = spattern.idx_position(0, 25);
  const std::tuple<Uint, Uint> res1 = spattern.idx_position(0, 31);

  BOOST_CHECK_EQUAL(std::get<0>(res0), 0u);
  BOOST_CHECK_EQUAL(std::get<1>(res0), 0u);

  BOOST_CHECK_EQUAL(std::get<0>(res1), 0u);
  BOOST_CHECK_EQUAL(std::get<1>(res1), 1u);

  // Check the second row
  const std::tuple<Uint, Uint> res2 = spattern.idx_position(41, 31);
  const std::tuple<Uint, Uint> res3 = spattern.idx_position(42, 110);

  BOOST_CHECK_EQUAL(std::get<0>(res2), 1u);
  BOOST_CHECK_EQUAL(std::get<1>(res2), 1u);

  BOOST_CHECK_EQUAL(std::get<0>(res3), 1u);
  BOOST_CHECK_EQUAL(std::get<1>(res3), 2u);

  // Check the third row
  const std::tuple<Uint, Uint> res4 = spattern.idx_position(107, 92);

  BOOST_CHECK_EQUAL(std::get<0>(res4), 2u);
  BOOST_CHECK_EQUAL(std::get<1>(res4), 2u);

  // spattern.print_vtu("block_sparsity.vtu");
}

// ----------------------------------------------------------------------------
