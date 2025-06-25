/// Generate automatically the 'main' function for the test module
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE block_multi_array_test
#include <boost/test/unit_test.hpp>

/// Standard template library headers
#include <cmath>
#include <ctime>
#include <iostream>

#include "common/BlockMultiArray.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

void build_block_array(common::BlockMultiArray<Real, Uint, Int> &block_array)
{
  std::unique_ptr<std::vector<Uint>> block_sizes(new std::vector<Uint>());
  (*block_sizes) = {2, 3, 4};

  std::unique_ptr<std::vector<Real>> field0_values(new std::vector<Real>());
  (*field0_values) = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};

  std::unique_ptr<std::vector<Uint>> field1_values(new std::vector<Uint>());
  (*field1_values) = {11, 12, 13, 14, 15, 16, 17, 18, 19};

  std::unique_ptr<std::vector<Int>> field2_values(new std::vector<Int>());
  (*field2_values) = {29, 28, 27, 26, 25, 24, 23, 22, 21};

  std::tuple<std::unique_ptr<std::vector<Real>>, std::unique_ptr<std::vector<Uint>>,
             std::unique_ptr<std::vector<Int>>>
      values;
  std::get<0>(values) = std::move(field0_values);
  std::get<1>(values) = std::move(field1_values);
  std::get<2>(values) = std::move(field2_values);

  block_array.build(std::move(values), std::move(block_sizes));
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_multi_array_construction)
{
  common::BlockMultiArray<Real, Uint, Int> block_array;
  BOOST_CHECK_EQUAL(block_array.empty(), true);

  build_block_array(block_array);

  common::BlockMultiArray<Real, Uint, Int> block_array2 = block_array;

  BOOST_CHECK_EQUAL(block_array2.size(), 9u);
  BOOST_CHECK_EQUAL(block_array2.nb_blocks(), 3u);
  BOOST_CHECK_EQUAL(block_array.empty(), false);

  // Check values in field 0
  const common::ArrayView<const Real, _1D, Uint> block0_0 = block_array2.const_block<0>(0);
  const common::ArrayView<const Real, _1D, Uint> block0_1 = block_array2.const_block<0>(1);
  const common::ArrayView<const Real, _1D, Uint> block0_2 = block_array2.const_block<0>(2);

  BOOST_CHECK_EQUAL(block0_0.size(), 2u);
  BOOST_CHECK_EQUAL(block0_1.size(), 3u);
  BOOST_CHECK_EQUAL(block0_2.size(), 4u);

  BOOST_CHECK_EQUAL(block0_0[0], 1.0);
  BOOST_CHECK_EQUAL(block0_0[1], 2.0);

  BOOST_CHECK_EQUAL(block0_1[0], 3.0);
  BOOST_CHECK_EQUAL(block0_1[1], 4.0);
  BOOST_CHECK_EQUAL(block0_1[2], 5.0);

  BOOST_CHECK_EQUAL(block0_2[0], 6.0);
  BOOST_CHECK_EQUAL(block0_2[1], 7.0);
  BOOST_CHECK_EQUAL(block0_2[2], 8.0);
  BOOST_CHECK_EQUAL(block0_2[3], 9.0);

  // Check values in field 1
  const common::ArrayView<const Uint, _1D, Uint> block1_0 = block_array2.const_block<1>(0);
  const common::ArrayView<const Uint, _1D, Uint> block1_1 = block_array2.const_block<1>(1);
  const common::ArrayView<const Uint, _1D, Uint> block1_2 = block_array2.const_block<1>(2);

  BOOST_CHECK_EQUAL(block1_0.size(), 2u);
  BOOST_CHECK_EQUAL(block1_1.size(), 3u);
  BOOST_CHECK_EQUAL(block1_2.size(), 4u);

  BOOST_CHECK_EQUAL(block1_0[0], 11u);
  BOOST_CHECK_EQUAL(block1_0[1], 12u);

  BOOST_CHECK_EQUAL(block1_1[0], 13u);
  BOOST_CHECK_EQUAL(block1_1[1], 14u);
  BOOST_CHECK_EQUAL(block1_1[2], 15u);

  BOOST_CHECK_EQUAL(block1_2[0], 16u);
  BOOST_CHECK_EQUAL(block1_2[1], 17u);
  BOOST_CHECK_EQUAL(block1_2[2], 18u);
  BOOST_CHECK_EQUAL(block1_2[3], 19u);

  // Check values in field 1
  const common::ArrayView<const Int, _1D, Uint> block2_0 = block_array2.const_block<2>(0);
  const common::ArrayView<const Int, _1D, Uint> block2_1 = block_array2.const_block<2>(1);
  const common::ArrayView<const Int, _1D, Uint> block2_2 = block_array2.const_block<2>(2);

  BOOST_CHECK_EQUAL(block2_0.size(), 2u);
  BOOST_CHECK_EQUAL(block2_1.size(), 3u);
  BOOST_CHECK_EQUAL(block2_2.size(), 4u);

  BOOST_CHECK_EQUAL(block2_0[0], 29);
  BOOST_CHECK_EQUAL(block2_0[1], 28);

  BOOST_CHECK_EQUAL(block2_1[0], 27);
  BOOST_CHECK_EQUAL(block2_1[1], 26);
  BOOST_CHECK_EQUAL(block2_1[2], 25);

  BOOST_CHECK_EQUAL(block2_2[0], 24);
  BOOST_CHECK_EQUAL(block2_2[1], 23);
  BOOST_CHECK_EQUAL(block2_2[2], 22);
  BOOST_CHECK_EQUAL(block2_2[3], 21);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_multi_array_functions)
{
  common::BlockMultiArray<Real, Uint, Int> block_array;
  build_block_array(block_array);

  const common::Range1D<Uint> limits0 = block_array.block_limits(0);
  BOOST_CHECK_EQUAL(limits0.lbound(), 0u);
  BOOST_CHECK_EQUAL(limits0.ubound(), 1u);

  const common::Range1D<Uint> limits1 = block_array.block_limits(1);
  BOOST_CHECK_EQUAL(limits1.lbound(), 2u);
  BOOST_CHECK_EQUAL(limits1.ubound(), 4u);

  const common::Range1D<Uint> limits2 = block_array.block_limits(2);
  BOOST_CHECK_EQUAL(limits2.lbound(), 5u);
  BOOST_CHECK_EQUAL(limits2.ubound(), 8u);

  // Check 'front_block' function
  common::ArrayView<const Real, _1D, Uint> front_block_f0 = block_array.front_block<0>();

  BOOST_CHECK_EQUAL(front_block_f0.size(), 2u);
  BOOST_CHECK_EQUAL(front_block_f0[0], 1.0);
  BOOST_CHECK_EQUAL(front_block_f0[1], 2.0);

  common::ArrayView<const Uint, _1D, Uint> front_block_f1 = block_array.front_block<1>();

  BOOST_CHECK_EQUAL(front_block_f1.size(), 2u);
  BOOST_CHECK_EQUAL(front_block_f1[0], 11u);
  BOOST_CHECK_EQUAL(front_block_f1[1], 12u);

  common::ArrayView<const Int, _1D, Uint> front_block_f2 = block_array.front_block<2>();

  BOOST_CHECK_EQUAL(front_block_f2.size(), 2u);
  BOOST_CHECK_EQUAL(front_block_f2[0], 29);
  BOOST_CHECK_EQUAL(front_block_f2[1], 28);

  // Check 'back_block' function
  common::ArrayView<const Real, _1D, Uint> back_block_f0 = block_array.back_block<0>();

  BOOST_CHECK_EQUAL(back_block_f0.size(), 4u);
  BOOST_CHECK_EQUAL(back_block_f0[0], 6.0);
  BOOST_CHECK_EQUAL(back_block_f0[1], 7.0);
  BOOST_CHECK_EQUAL(back_block_f0[2], 8.0);
  BOOST_CHECK_EQUAL(back_block_f0[3], 9.0);

  common::ArrayView<const Uint, _1D, Uint> back_block_f1 = block_array.back_block<1>();

  BOOST_CHECK_EQUAL(back_block_f1.size(), 4u);
  BOOST_CHECK_EQUAL(back_block_f1[0], 16u);
  BOOST_CHECK_EQUAL(back_block_f1[1], 17u);
  BOOST_CHECK_EQUAL(back_block_f1[2], 18u);
  BOOST_CHECK_EQUAL(back_block_f1[3], 19u);

  common::ArrayView<const Int, _1D, Uint> back_block_f2 = block_array.back_block<2>();

  BOOST_CHECK_EQUAL(back_block_f2.size(), 4u);
  BOOST_CHECK_EQUAL(back_block_f2[0], 24);
  BOOST_CHECK_EQUAL(back_block_f2[1], 23);
  BOOST_CHECK_EQUAL(back_block_f2[2], 22);
  BOOST_CHECK_EQUAL(back_block_f2[3], 21);

  common::ArrayView<Int, _1D, Uint> mid_block_f2 = block_array.block<2>(1);
  BOOST_CHECK_EQUAL(mid_block_f2.size(), 3u);
  mid_block_f2[0] = 100;
  mid_block_f2[1] = 101;
  mid_block_f2[2] = 102;

  const common::ArrayView<const Int, _1D, Uint> const_mid_block_f2 = block_array.const_block<2>(1);
  BOOST_CHECK_EQUAL(const_mid_block_f2.size(), 3u);
  BOOST_CHECK_EQUAL(const_mid_block_f2[0], 100);
  BOOST_CHECK_EQUAL(const_mid_block_f2[1], 101);
  BOOST_CHECK_EQUAL(const_mid_block_f2[2], 102);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_multi_array_block_sorting)
{
  common::BlockMultiArray<Real, Uint, Int> block_array;
  build_block_array(block_array);

  block_array.sort_blocks<2>();

  // --------------

  common::ArrayView<const Real, _1D, Uint> b0_f0 = block_array.const_block<0>(0);
  BOOST_CHECK_EQUAL(b0_f0.size(), 2u);
  BOOST_CHECK_EQUAL(b0_f0[0], 2.0);
  BOOST_CHECK_EQUAL(b0_f0[1], 1.0);

  common::ArrayView<const Real, _1D, Uint> b1_f0 = block_array.const_block<0>(1);
  BOOST_CHECK_EQUAL(b1_f0.size(), 3u);
  BOOST_CHECK_EQUAL(b1_f0[0], 5.0);
  BOOST_CHECK_EQUAL(b1_f0[1], 4.0);
  BOOST_CHECK_EQUAL(b1_f0[2], 3.0);

  common::ArrayView<const Real, _1D, Uint> b2_f0 = block_array.const_block<0>(2);
  BOOST_CHECK_EQUAL(b2_f0.size(), 4u);
  BOOST_CHECK_EQUAL(b2_f0[0], 9.0);
  BOOST_CHECK_EQUAL(b2_f0[1], 8.0);
  BOOST_CHECK_EQUAL(b2_f0[2], 7.0);
  BOOST_CHECK_EQUAL(b2_f0[3], 6.0);

  // --------------

  common::ArrayView<const Uint, _1D, Uint> b0_f1 = block_array.const_block<1>(0);
  BOOST_CHECK_EQUAL(b0_f1.size(), 2u);
  BOOST_CHECK_EQUAL(b0_f1[0], 12);
  BOOST_CHECK_EQUAL(b0_f1[1], 11);

  common::ArrayView<const Uint, _1D, Uint> b1_f1 = block_array.const_block<1>(1);
  BOOST_CHECK_EQUAL(b1_f1.size(), 3u);
  BOOST_CHECK_EQUAL(b1_f1[0], 15);
  BOOST_CHECK_EQUAL(b1_f1[1], 14);
  BOOST_CHECK_EQUAL(b1_f1[2], 13);

  common::ArrayView<const Uint, _1D, Uint> b2_f1 = block_array.const_block<1>(2);
  BOOST_CHECK_EQUAL(b2_f1.size(), 4u);
  BOOST_CHECK_EQUAL(b2_f1[0], 19);
  BOOST_CHECK_EQUAL(b2_f1[1], 18);
  BOOST_CHECK_EQUAL(b2_f1[2], 17);
  BOOST_CHECK_EQUAL(b2_f1[3], 16);

  // --------------

  common::ArrayView<const Int, _1D, Uint> b0_f2 = block_array.const_block<2>(0);
  BOOST_CHECK_EQUAL(b0_f2.size(), 2u);
  BOOST_CHECK_EQUAL(b0_f2[0], 28);
  BOOST_CHECK_EQUAL(b0_f2[1], 29);

  common::ArrayView<const Int, _1D, Uint> b1_f2 = block_array.const_block<2>(1);
  BOOST_CHECK_EQUAL(b1_f2.size(), 3u);
  BOOST_CHECK_EQUAL(b1_f2[0], 25);
  BOOST_CHECK_EQUAL(b1_f2[1], 26);
  BOOST_CHECK_EQUAL(b1_f2[2], 27);

  common::ArrayView<const Int, _1D, Uint> b2_f2 = block_array.const_block<2>(2);
  BOOST_CHECK_EQUAL(b2_f2.size(), 4u);
  BOOST_CHECK_EQUAL(b2_f2[0], 21);
  BOOST_CHECK_EQUAL(b2_f2[1], 22);
  BOOST_CHECK_EQUAL(b2_f2[2], 23);
  BOOST_CHECK_EQUAL(b2_f2[3], 24);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_multi_array_block_resize)
{
  common::BlockMultiArray<Real, Uint, Int> block_array;
  build_block_array(block_array);

  block_array.resize_blocks({1, 3, 6});

  BOOST_CHECK_EQUAL(block_array.size(), 10u);
  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 3u);

  const common::ArrayView<const Real, _1D, Uint> block0 = block_array.const_block<0>(0);
  const common::ArrayView<const Real, _1D, Uint> block1 = block_array.const_block<0>(1);
  const common::ArrayView<const Real, _1D, Uint> block2 = block_array.const_block<0>(2);

  BOOST_CHECK_EQUAL(block0.size(), 1u);
  BOOST_CHECK_EQUAL(block1.size(), 3u);
  BOOST_CHECK_EQUAL(block2.size(), 6u);

  common::ArrayView<const Real, _1D, Uint> back_block_f0 = block_array.back_block<0>();

  BOOST_CHECK_EQUAL(back_block_f0.size(), 6u);
  BOOST_CHECK_EQUAL(back_block_f0[0], 6.0);
  BOOST_CHECK_EQUAL(back_block_f0[1], 7.0);
  BOOST_CHECK_EQUAL(back_block_f0[2], 8.0);
  BOOST_CHECK_EQUAL(back_block_f0[3], 9.0);
  BOOST_CHECK_EQUAL(back_block_f0[4], 0.0);
  BOOST_CHECK_EQUAL(back_block_f0[5], 0.0);

  common::ArrayView<const Uint, _1D, Uint> back_block_f1 = block_array.back_block<1>();

  BOOST_CHECK_EQUAL(back_block_f1.size(), 6u);
  BOOST_CHECK_EQUAL(back_block_f1[0], 16u);
  BOOST_CHECK_EQUAL(back_block_f1[1], 17u);
  BOOST_CHECK_EQUAL(back_block_f1[2], 18u);
  BOOST_CHECK_EQUAL(back_block_f1[3], 19u);
  BOOST_CHECK_EQUAL(back_block_f1[4], 0u);
  BOOST_CHECK_EQUAL(back_block_f1[5], 0u);

  common::ArrayView<const Int, _1D, Uint> back_block_f2 = block_array.back_block<2>();

  BOOST_CHECK_EQUAL(back_block_f2.size(), 6u);
  BOOST_CHECK_EQUAL(back_block_f2[0], 24);
  BOOST_CHECK_EQUAL(back_block_f2[1], 23);
  BOOST_CHECK_EQUAL(back_block_f2[2], 22);
  BOOST_CHECK_EQUAL(back_block_f2[3], 21);
  BOOST_CHECK_EQUAL(back_block_f2[4], 0);
  BOOST_CHECK_EQUAL(back_block_f2[5], 0);
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(block_multi_array_block_remove)
{
  common::BlockMultiArray<Real, Uint, Int> block_array;
  build_block_array(block_array);

  // Remove the first and third blocks, leaving only the middle one ...
  block_array.remove_blocks({0, 2});

  BOOST_CHECK_EQUAL(block_array.size(), 3u);
  BOOST_CHECK_EQUAL(block_array.nb_blocks(), 1u);

  const common::ArrayView<const Real, _1D, Uint> block_f0 = block_array.const_block<0>(0);
  const common::ArrayView<const Uint, _1D, Uint> block_f1 = block_array.const_block<1>(0);
  const common::ArrayView<const Int, _1D, Uint> block_f2  = block_array.const_block<2>(0);

  BOOST_CHECK_EQUAL(block_f0.size(), 3u);
  BOOST_CHECK_EQUAL(block_f0[0], 3.0);
  BOOST_CHECK_EQUAL(block_f0[1], 4.0);
  BOOST_CHECK_EQUAL(block_f0[2], 5.0);

  BOOST_CHECK_EQUAL(block_f1.size(), 3u);
  BOOST_CHECK_EQUAL(block_f1[0], 13.0);
  BOOST_CHECK_EQUAL(block_f1[1], 14.0);
  BOOST_CHECK_EQUAL(block_f1[2], 15.0);

  BOOST_CHECK_EQUAL(block_f2.size(), 3u);
  BOOST_CHECK_EQUAL(block_f2[0], 27.0);
  BOOST_CHECK_EQUAL(block_f2[1], 26.0);
  BOOST_CHECK_EQUAL(block_f2[2], 25.0);
}

// ----------------------------------------------------------------------------
