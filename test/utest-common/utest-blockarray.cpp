/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE block_array_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

#include "common/BlockArray.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_array_construction)
{
  common::BlockArray<Real, Uint> block_array;
  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(block_array.size(), 0u);

  std::vector<Real> raw_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  block_array.reserve(9, 3);
  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data(), 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 3, 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 6, 3));

  std::cout << block_array << std::endl;

  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 3u);
  BOOST_CHECK_EQUAL(block_array.size(), 9u);

  common::BlockArray<Real, Uint> block_array2;

  std::unique_ptr<std::vector<Real>> values(new std::vector<Real>());
  std::unique_ptr<std::vector<Uint>> block_sizes(new std::vector<Uint>());

  (*values)      = {101, 102, 103, 104, 105, 106, 107, 108, 109};
  (*block_sizes) = {2, 4, 3};

  block_array2.build(std::move(values), std::move(block_sizes));

  BOOST_CHECK_EQUAL(block_array2.nb_blocks(), 3u);
  BOOST_CHECK_EQUAL(block_array2.size(), 9u);

  std::cout << block_array2 << std::endl;

  std::vector<Uint> new_block_sizes = {3, 3, 7};
  block_array2.resize_blocks(new_block_sizes);
  const common::ArrayView<const Real, _1D, Uint> block0 = block_array2.const_block(0);
  const common::ArrayView<const Real, _1D, Uint> block1 = block_array2.const_block(1);
  const common::ArrayView<const Real, _1D, Uint> block2 = block_array2.const_block(2);

  BOOST_CHECK_EQUAL(block_array2.size(), 13u);

  BOOST_CHECK_EQUAL(block0.size(), 3u);
  BOOST_CHECK_EQUAL(block1.size(), 3u);
  BOOST_CHECK_EQUAL(block2.size(), 7u);

  BOOST_CHECK_EQUAL(block0[0], 101.0);
  BOOST_CHECK_EQUAL(block0[1], 102.0);
  BOOST_CHECK_EQUAL(block0[2], 0.0);

  BOOST_CHECK_EQUAL(block1[0], 103.0);
  BOOST_CHECK_EQUAL(block1[1], 104.0);
  BOOST_CHECK_EQUAL(block1[2], 105.0);

  BOOST_CHECK_EQUAL(block2[0], 107.0);
  BOOST_CHECK_EQUAL(block2[1], 108.0);
  BOOST_CHECK_EQUAL(block2[2], 109.0);
  BOOST_CHECK_EQUAL(block2[3], 0.0);
  BOOST_CHECK_EQUAL(block2[4], 0.0);
  BOOST_CHECK_EQUAL(block2[5], 0.0);
  BOOST_CHECK_EQUAL(block2[6], 0.0);

  std::cout << block_array2 << std::endl;

  block_array2.remove_blocks({1});

  BOOST_CHECK_EQUAL(block_array2.nb_blocks(), 2u);
  BOOST_CHECK_EQUAL(block_array2.size(), 10u);

  const common::ArrayView<const Real, _1D, Uint> blk0 = block_array2.const_block(0);
  const common::ArrayView<const Real, _1D, Uint> blk1 = block_array2.const_block(1);

  BOOST_CHECK_EQUAL(blk0.size(), 3u);
  BOOST_CHECK_EQUAL(blk1.size(), 7u);

  BOOST_CHECK_EQUAL(blk0[0], 101.0);
  BOOST_CHECK_EQUAL(blk0[1], 102.0);
  BOOST_CHECK_EQUAL(blk0[2], 0.0);

  BOOST_CHECK_EQUAL(blk1[0], 107.0);
  BOOST_CHECK_EQUAL(blk1[1], 108.0);
  BOOST_CHECK_EQUAL(blk1[2], 109.0);
  BOOST_CHECK_EQUAL(blk1[3], 0.0);
  BOOST_CHECK_EQUAL(blk1[4], 0.0);
  BOOST_CHECK_EQUAL(blk1[5], 0.0);
  BOOST_CHECK_EQUAL(blk1[6], 0.0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_array_io)
{
  common::BlockArray<Real, Uint> block_array;
  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(block_array.size(), 0u);

  std::vector<Real> raw_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  block_array.reserve(9, 3);

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data(), 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 3, 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 6, 3));

  auto b1 = block_array.const_block(1);
  std::cout << "Block size = " << b1.size() << std::endl;
  std::cout << "Block values = " << std::endl;
  for (Uint i = 0; i < b1.size(); ++i)
  {
    std::cout << b1[i] << " ";
  }
  std::cout << std::endl;

  std::cout << "Last block:" << std::endl;
  auto last_block = block_array.back_block();

  std::cout << "Block size = " << last_block.size() << std::endl;
  std::cout << "Block values = " << std::endl;
  for (Uint i = 0; i < last_block.size(); ++i)
  {
    std::cout << last_block[i] << " ";
  }
  std::cout << std::endl;

  raw_values[3] = 20;
  raw_values[4] = 40;
  raw_values[5] = 30;
  block_array.insert_block(2, common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 3, 3));

  std::cout << "Block array before sort: " << std::endl;
  std::cout << block_array << std::endl;

  struct Compare
  {
    bool operator()(const Real a, const Real b)
    {
      return a > b;
    }
  };

  Compare comp;

  block_array.sort_blocks(comp);
  // block_array.sort_blocks<Compare>();
  // block_array.sort_blocks();
  std::cout << "Block array after sort: " << std::endl;
  std::cout << block_array << std::endl;
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_array_raw_data_access)
{
  common::BlockArray<Real, Uint> block_array;
  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 0u);
  BOOST_CHECK_EQUAL(block_array.size(), 0u);

  std::vector<Real> raw_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  block_array.reserve(9, 3);

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data(), 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 3, 3));

  block_array.create_back_block(3);
  block_array.fill_last_block(common::ArrayView<const Real, _1D, Uint>(raw_values.data() + 6, 3));

  common::VectorPack<const std::vector<Real>, const std::vector<Uint>> raw_crefs =
      block_array.raw_data_crefs();

  const auto &raw_values_ref  = raw_crefs.get<0>().get();
  const auto &raw_offsets_ref = raw_crefs.get<1>().get();

  std::cout << raw_values_ref << std::endl;
  std::cout << raw_offsets_ref << std::endl;

  BOOST_CHECK_EQUAL(raw_values_ref.size(), 9);
  for (Uint i = 0; i < raw_values_ref.size(); ++i)
  {
    BOOST_CHECK_EQUAL(raw_values_ref[i], raw_values[i]);
  }

  BOOST_CHECK_EQUAL(raw_offsets_ref.size(), 4);
  BOOST_CHECK_EQUAL(raw_offsets_ref[0], 0);
  BOOST_CHECK_EQUAL(raw_offsets_ref[1], 3);
  BOOST_CHECK_EQUAL(raw_offsets_ref[2], 6);
  BOOST_CHECK_EQUAL(raw_offsets_ref[3], 9);
}

// ----------------------------------------------------------------------------
