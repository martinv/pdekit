/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE tags_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iostream>

/// PDEKIT headers
#include "common/Constants.hpp"
#include "common/Meta.hpp"

using namespace pdekit;
using namespace boost::unit_test;

BOOST_AUTO_TEST_CASE(tag_data_structures_utest)
{

  // Test that some metaprogramming utilities used to define the static
  // variables in the tag class are correct
  BOOST_CHECK_EQUAL(common::PowerOfTwo<1>::value, 2u);       // 2^1 = 2
  BOOST_CHECK_EQUAL(common::PowerOfTwo<2>::value, 4u);       // 2^2 = 4
  BOOST_CHECK_EQUAL(common::PowerOfTwo<3>::value, 8u);       // 2^3 = 8
  BOOST_CHECK_EQUAL(common::PowerOfTwo<4>::value, 16u);      // 2^4 = 16
  BOOST_CHECK_EQUAL(common::PowerOfTwo<5>::value, 32u);      // 2^5 = 32
  BOOST_CHECK_EQUAL(common::PowerOfTwo<6>::value, 64u);      // 2^6 = 64
  BOOST_CHECK_EQUAL(common::PowerOfTwo<7>::value, 128u);     // 2^7 = 128
  BOOST_CHECK_EQUAL(common::PowerOfTwo<8>::value, 256u);     // 2^8 = 256
  BOOST_CHECK_EQUAL(common::PowerOfTwo<9>::value, 512u);     // 2^9 = 512
  BOOST_CHECK_EQUAL(common::PowerOfTwo<10>::value, 1024u);   // 2^10 = 1024
  BOOST_CHECK_EQUAL(common::PowerOfTwo<11>::value, 2048u);   // 2^11 = 2048
  BOOST_CHECK_EQUAL(common::PowerOfTwo<12>::value, 4096u);   // 2^12 = 4096
  BOOST_CHECK_EQUAL(common::PowerOfTwo<13>::value, 8192u);   // 2^13 = 8192
  BOOST_CHECK_EQUAL(common::PowerOfTwo<14>::value, 16384u);  // 2^14 = 16384
  BOOST_CHECK_EQUAL(common::PowerOfTwo<15>::value, 32768u);  // 2^15 = 32768
  BOOST_CHECK_EQUAL(common::PowerOfTwo<16>::value, 65536u);  // 2^16 = 65536
  BOOST_CHECK_EQUAL(common::PowerOfTwo<17>::value, 131072u); // 2^17 = 131072
  BOOST_CHECK_EQUAL(common::PowerOfTwo<18>::value, 262144u); // 2^18 = 261144

  BOOST_CHECK_EQUAL(common::MinNbBitsToStoreNumber<2>::value, 2);
  BOOST_CHECK_EQUAL(common::MinNbBitsToStoreNumber<3>::value, 2);
  BOOST_CHECK_EQUAL(common::MinNbBitsToStoreNumber<4>::value, 3);
  BOOST_CHECK_EQUAL(common::MinNbBitsToStoreNumber<5>::value, 3);
  BOOST_CHECK_EQUAL(common::MinNbBitsToStoreNumber<15>::value, 4);
}

/*
BOOST_AUTO_TEST_CASE(tag_test_case)
{
  common::Tag3<ElemShape, PolyOrder, RefTopology> tag_a, tag_b, tag_c;

  tag_a.set_key<ElemShape>(Triag);
  tag_a.set_key<PolyOrder>(P4);
  tag_a.set_key<RefTopology>(Equidist);
  BOOST_CHECK_EQUAL(tag_a.key<ElemShape>(), Triag);
  BOOST_CHECK_EQUAL(tag_a.key<PolyOrder>(), P4);
  BOOST_CHECK_EQUAL(tag_a.key<RefTopology>(), Equidist);

  std::cout << common::Tag3<ElemShape, PolyOrder,
RefTopology>::keys_to_string(Triag, P4, Equidist)
            << std::endl;

  tag_b = common::Tag3<ElemShape, PolyOrder,
RefTopology>::string_to_tag("Triag-P4-Equidist");

  BOOST_CHECK_EQUAL(tag_a.as_string(), tag_b.as_string());

  std::cout << tag_b.as_string() << std::endl;

  tag_c = common::Tag3<ElemShape, PolyOrder, RefTopology>::make_tag(Hexa, P6,
Warpblend);

  BOOST_CHECK_EQUAL(tag_c.key<ElemShape>(), Hexa);
  BOOST_CHECK_EQUAL(tag_c.key<PolyOrder>(), P6);
  BOOST_CHECK_EQUAL(tag_c.key<RefTopology>(), Warpblend);

  tag_a = common::Tag3<ElemShape, PolyOrder, RefTopology>(
      tag_c.key<ElemShape>(), tag_c.key<PolyOrder>(), tag_c.key<RefTopology>());

  BOOST_CHECK_EQUAL(tag_a.key<ElemShape>(), tag_c.key<ElemShape>());
  BOOST_CHECK_EQUAL(tag_a.key<PolyOrder>(), tag_c.key<PolyOrder>());
  BOOST_CHECK_EQUAL(tag_a.key<RefTopology>(), tag_c.key<RefTopology>());

  BOOST_CHECK_EQUAL(tag_a == tag_c, true);
}
*/
