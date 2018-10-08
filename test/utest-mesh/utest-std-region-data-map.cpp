/// Generate automatically the 'main' function for the test module
#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE data_map_utest
#include <boost/test/unit_test.hpp>

/// STL headers
#include <iomanip>
#include <iostream>

#include "common/DataMap.hpp"
#include "interpolation/FEValues.hpp"
#include "mesh/std_region/PointSetTagExt.hpp"

using namespace pdekit;

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(check_std_region_data_map_utest)
{

  common::DataMap<mesh::PointSetTagExt, math::DenseDMat<Real>> map1;

  BOOST_CHECK_EQUAL(map1.size(), 0U);

  for (Uint P = 1; P <= 5; ++P)
  {
    map1.create(mesh::PointSetTagExt(mesh::PointSetTag(ElemShape::Quad, P, PointSetID::Equidist),
                                     P0, mesh::CellTransform::NO_TRANS, 0));

    const mesh::PointSetTag std_reg_tag =
        mesh::PointSetTag(ElemShape::Quad, P, PointSetID::Equidist);

    common::PtrHandle<math::DenseDMat<Real>> entry = map1.std_region_data(
        mesh::PointSetTagExt(std_reg_tag, P0, mesh::CellTransform::NO_TRANS, 0));

    (*entry).resize(10, 10);
  }

  common::PtrHandle<math::DenseDMat<Real>> entry = map1.std_region_data(
      mesh::PointSetTagExt(mesh::PointSetTag(ElemShape::Quad, P1, PointSetID::Equidist), P0,
                           mesh::CellTransform::NO_TRANS, 0));

  for (common::DataMap<mesh::PointSetTagExt, math::DenseDMat<Real>>::const_iterator it =
           map1.cbegin();
       it != map1.cend(); ++it)
  {
    std::cout << it.key_value().std_region_tag().as_string() << std::endl;
    std::cout << (*it.data_ptr());
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(custom_key_std_data_cache_utest)
{

  typedef std::pair<mesh::PointSetTag, Uint> CustomKey;

  struct CustomHasher
  {
    std::size_t operator()(const CustomKey &key) const
    {
      return key.first.store_value() ^ key.second;
    }
  };

  struct CustomDataMapTrait
  {
    using hash_type   = CustomHasher;
    using container_t = common::DataMapStorageByMap;
  };

  common::DataMap<CustomKey, math::DenseDMat<Real>, CustomDataMapTrait> custom_map;

  BOOST_CHECK_EQUAL(custom_map.size(), 0U);

  for (Uint P = 1; P <= 5; ++P)
  {
    custom_map.create(
        std::make_pair(mesh::PointSetTag(ElemShape::Triag, P, PointSetID::Equidist), 0u));

    const mesh::PointSetTag std_reg_tag =
        mesh::PointSetTag(ElemShape::Triag, P, PointSetID::Equidist);

    common::PtrHandle<math::DenseDMat<Real>> entry =
        custom_map.std_region_data(std::pair<mesh::PointSetTag, Uint>(std_reg_tag, 0u));

    (*entry).resize(10, 10);
  }

  for (common::DataMap<CustomKey, math::DenseDMat<Real>, CustomDataMapTrait>::const_iterator it =
           custom_map.cbegin();
       it != custom_map.cend(); ++it)
  {
    std::cout << it.key_value().first.as_string() << std::endl;
    std::cout << (*it.data_ptr());
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(key_pair_std_data_cache_utest)
{
  typedef std::pair<math::DenseDMat<Real>, math::DenseDMat<Real>> value_type;
  typedef mesh::PointSetTagExtPair key_type;

  struct CustomDataMapTrait
  {
    using hash_type   = mesh::detail::PointSetTagExtPairHasher;
    using container_t = common::DataMapStorageByMap;
  };

  common::DataMap<key_type, value_type, CustomDataMapTrait> pair_map;

  for (Uint P = 1; P <= 5; ++P)
  {
    const mesh::PointSetTag t1(ElemShape::Triag, P, PointSetID::Equidist);
    const mesh::PointSetTag t2(ElemShape::Triag, P + 1, PointSetID::Equidist);

    key_type key(
        mesh::PointSetTagExtPair(mesh::PointSetTagExt(t1, 0, mesh::CellTransform::NO_TRANS, 0),
                                 mesh::PointSetTagExt(t2, 0, mesh::CellTransform::NO_TRANS, 0)));

    pair_map.create(key);

    common::PtrHandle<std::pair<math::DenseDMat<Real>, math::DenseDMat<Real>>> entry =
        pair_map.std_region_data(key);

    (*entry).first.resize(5, 5);
    (*entry).second.resize(10, 10);
  }

  BOOST_CHECK_EQUAL(pair_map.size(), 5u);

  for (common::DataMap<key_type, value_type, CustomDataMapTrait>::const_iterator it =
           pair_map.cbegin();
       it != pair_map.cend(); ++it)
  {
    std::cout << it.key_value().m_t1.as_string() << " : " << it.key_value().m_t2.as_string()
              << std::endl;
    BOOST_CHECK_EQUAL((*it.data_ptr()).first.rows(), 5u);
    BOOST_CHECK_EQUAL((*it.data_ptr()).first.cols(), 5u);
    BOOST_CHECK_EQUAL((*it.data_ptr()).second.rows(), 10u);
    BOOST_CHECK_EQUAL((*it.data_ptr()).second.cols(), 10u);

    std::cout << (*it.data_ptr()).first << std::endl;
    std::cout << (*it.data_ptr()).second;
    std::cout << std::endl;
  }
}

// ----------------------------------------------------------------------------
