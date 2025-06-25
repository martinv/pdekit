/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE vector_utils_utest
#include <boost/bind.hpp>
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/VectorUtils.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vector_print_utest)
{
  std::vector<Real> v1 = {1.1, 2.1, 3.1};
  std::vector<Int> v2  = {10, 11, 12, 13};

  std::cout << "[" << v1 << "]" << std::endl;
  std::cout << "[" << v2 << "]" << std::endl;
}
// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(vector_pack_utest)
{
  std::vector<Real> v1 = {1.1, 2.1, 3.1};
  std::vector<Int> v2  = {10, 11, 12, 13};

  // VectorPack<std::vector<Real>, std::vector<Int>> vp(v1, v2);
  auto vp = common::pack_vectors(v1, v2);

  // auto &v1_ref = vp.get<0>().get();
  common::VectorPackEntryType<0, std::vector<Real>, std::vector<Int>>::type &v1_ref =
      vp.get<0>().get();
  auto &v2_ref = vp.get<1>().get();

  v1_ref[0] = 100.0;
  v1_ref[1] = 201.0;
  v1_ref[2] = 300.0;

  v2_ref[2] = 3;

  BOOST_CHECK_EQUAL(v1_ref.size(), 3);
  BOOST_CHECK_EQUAL(v2_ref.size(), 4);

  BOOST_CHECK_EQUAL(v1[0], 100.0);
  BOOST_CHECK_EQUAL(v1[1], 201.0);
  BOOST_CHECK_EQUAL(v1[2], 300.0);

  std::cout << "[" << v1 << "]" << std::endl;
  std::cout << "[" << v2 << "]" << std::endl;
}

// ----------------------------------------------------------------------------
