/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE hash_utest
#include <boost/bind.hpp>
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/Hash.hpp"
#include "common/MPI/MPIEnv.hpp"

using namespace pdekit;
using namespace boost::unit_test;

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(hash_utest)
{
  std::vector<size_t> vec1 = {4, 5, 6, 7};

  const size_t hash_val1 = common::hash_seq_accum(vec1);
  std::cout << "Hashed value 1: " << hash_val1 << std::endl;

  std::vector<size_t> vec2 = {6, 7, 4, 5};

  const size_t hash_val2 = common::hash_seq_accum(vec2);
  std::cout << "Hashed value 2: " << hash_val2 << std::endl;
}

// ----------------------------------------------------------------------------
