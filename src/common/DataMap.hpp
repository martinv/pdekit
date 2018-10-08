#ifndef PDEKIT_Common_Data_Map_hpp
#define PDEKIT_Common_Data_Map_hpp

#include "common/DataMapImpl.hpp"
#include "common/Meta.hpp"

namespace pdekit
{

namespace common
{

struct DataMapStorageByVector
{
};
struct DataMapStorageByMap
{
};

struct DefaultDataMapTrait
{
  using container_t = DataMapStorageByVector;
  using hash_type   = NullType;
};

template <typename KeyType, typename T, typename Traits = DefaultDataMapTrait>
class DataMap
{
  private:
  /*
  // using map_implementation = detail::DataMapImplWithThreads<KeyType, T,
  Hasher>; using map_implementation = detail::VecDataMapImplNoThreads<KeyType,
  T>;
  */

  using map_implementation = typename SelectType<
      TypesAreIdentical<typename Traits::container_t, DataMapStorageByVector>::value,
      detail::VecDataMapImplNoThreads<KeyType, T>,
      typename SelectType<
          TypesAreIdentical<typename Traits::container_t, DataMapStorageByMap>::value,
          detail::DataMapImplWithThreads<KeyType, T, typename Traits::hash_type>,
          FalseType>::type>::type;

  public:
  /// TYPEDEFS

  using key_type       = typename map_implementation::key_type;
  using const_iterator = typename map_implementation::const_iterator;
  using iterator       = typename map_implementation::iterator;

  /// METHODS

  /// Default constructor
  DataMap();

  /// Default destructor
  ~DataMap();

  /// Return size of the underlying map
  Uint size() const;

  /// Remove all internal data
  void clear();

  /// Create a new entry
  common::PtrHandle<T> create(KeyType const key);

  /// Create a new entry, const version
  const common::PtrHandle<T const> create(KeyType const key) const;

  /// Get one entry
  common::PtrHandle<T> std_region_data(KeyType const key);

  /// Get one entry, const version
  const common::PtrHandle<T const> std_region_data(KeyType const key) const;

  /// Return the constant 'begin' iterator over the data
  inline const_iterator cbegin() const
  {
    return m_impl.cbegin();
  }

  /// Return the constant 'end' iterator over the data
  inline const_iterator cend() const
  {
    return m_impl.cend();
  }

  inline iterator begin()
  {
    return m_impl.begin();
  }

  inline iterator end()
  {
    return m_impl.end();
  }

  private:
  // The actual map implementation
  map_implementation m_impl;
};

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
DataMap<KeyType, T, Traits>::DataMap()
{
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
DataMap<KeyType, T, Traits>::~DataMap()
{
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
Uint DataMap<KeyType, T, Traits>::size() const
{
  return m_impl.size();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
void DataMap<KeyType, T, Traits>::clear()
{
  m_impl.clear();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
common::PtrHandle<T> DataMap<KeyType, T, Traits>::create(KeyType const key)
{
  return m_impl.create(key);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
const common::PtrHandle<T const> DataMap<KeyType, T, Traits>::create(KeyType const key) const
{
  return m_impl.create(key);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
common::PtrHandle<T> DataMap<KeyType, T, Traits>::std_region_data(KeyType const key)
{
  return m_impl.std_region_data(key);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Traits>
const common::PtrHandle<T const> DataMap<KeyType, T, Traits>::std_region_data(
    KeyType const key) const
{
  return m_impl.std_region_data(key);
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
