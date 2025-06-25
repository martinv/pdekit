#ifndef PDEKIT_Common_Data_Map_Impl_hpp
#define PDEKIT_Common_Data_Map_Impl_hpp

#include <mutex>
#include <unordered_map>
#include <vector>

#include "common/DataMapIterator.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------
// Implementation that does not rely on threads
// ----------------------------------------------------------------------------

namespace detail
{

template <typename KeyType, typename T, typename Hasher = typename KeyType::hash_type>
class DataMapImplNoThreads
{
  private:
  using const_iterator_impl = typename std::unordered_map<KeyType, T *, Hasher>::const_iterator;
  using iterator_impl       = typename std::unordered_map<KeyType, T *, Hasher>::iterator;

  public:
  /// TYPEDEFS

  using key_type = KeyType;

  using const_iterator = DataMapIterator<T, DataMapConstItTraits<T, KeyType, const_iterator_impl>>;

  using iterator = DataMapIterator<T, DataMapItTraits<T, KeyType, iterator_impl>>;

  /// METHODS

  /// Default constructor
  DataMapImplNoThreads();

  /// Default destructor
  ~DataMapImplNoThreads();

  /// Return size of the underlying map
  Uint size() const;

  /// Remove one entry
  // void remove(KeyType const key);

  /// Remove all internal data
  void clear();

  /// Create a new entry
  PtrHandle<T> create(KeyType const key);

  /// Create a new entry, const version
  const PtrHandle<T const> create(KeyType const key) const;

  /// Get one entry
  PtrHandle<T> std_region_data(KeyType const key);

  /// Get one entry, const version
  const PtrHandle<T const> std_region_data(KeyType const key) const;

  /// Return the constant 'begin' iterator over the data
  inline const_iterator cbegin() const
  {
    const_iterator it(m_std_region_data.cbegin());
    return it;
  }

  /// Return the constant 'end' iterator over the data
  inline const_iterator cend() const
  {
    const_iterator it(m_std_region_data.cend());
    return it;
  }

  inline iterator begin()
  {
    iterator it(m_std_region_data.begin());
    return it;
  }

  inline iterator end()
  {
    iterator it(m_std_region_data.end());
    return it;
  }

  private:
  /// TYPES
  /// To use a custom type as key in unordered map, two things have to be
  /// defined: 1) An equality operator:
  ///    const bool operator==(const Key &other) const
  ///    {
  ///      ...
  ///    }
  /// 2) A hash function (see DefaultStdRegMapKeyHasher at the top of this
  /// file)

  using raw_iterator       = typename std::unordered_map<KeyType, T *, Hasher>::iterator;
  using const_raw_iterator = typename std::unordered_map<KeyType, T *, Hasher>::const_iterator;

  /// Map [std. region id => data for this type of std region]
  std::unordered_map<KeyType, T *, Hasher> m_std_region_data;
  const_raw_iterator m_cached_entry;
};

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
DataMapImplNoThreads<KeyType, T, Hasher>::DataMapImplNoThreads()
{
  m_cached_entry = m_std_region_data.end();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
DataMapImplNoThreads<KeyType, T, Hasher>::~DataMapImplNoThreads()
{
  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete it->second;
    it->second = nullptr;
  }
  m_std_region_data.clear();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
Uint DataMapImplNoThreads<KeyType, T, Hasher>::size() const
{
  return m_std_region_data.size();
}

// ----------------------------------------------------------------------------

/*
template <typename KeyType, typename T, typename Hasher>
void StdRegionDataMapImplNoThreads<KeyType, T, Hasher>::remove(KeyType const
key)
{
  raw_iterator it = m_std_region_data.find(key);
  if (it != m_std_region_data.end())
  {
    delete it->second;
    m_std_region_data.erase(it);
  }

  if (!m_std_region_data.empty())
  {
    m_cached_entry = m_std_region_data.cbegin();
  }
  else
  {
    m_cached_entry = m_std_region_data.cend();
  }
}
*/

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
void DataMapImplNoThreads<KeyType, T, Hasher>::clear()
{
  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete it->second;
  }
  m_std_region_data.clear();
  m_cached_entry = m_std_region_data.cend();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
PtrHandle<T> DataMapImplNoThreads<KeyType, T, Hasher>::create(KeyType const key)
{
  raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T>(it->second);
  }

  T *new_data = new T();
  m_std_region_data.insert(std::pair<KeyType, T *>(key, new_data));
  return PtrHandle<T>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
const PtrHandle<T const> DataMapImplNoThreads<KeyType, T, Hasher>::create(KeyType const key) const
{
  const_raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T const>(it->second);
  }

  T *new_data = new T();
  m_std_region_data.insert(std::pair<KeyType, T *>(key, new_data));
  return PtrHandle<T const>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
PtrHandle<T> DataMapImplNoThreads<KeyType, T, Hasher>::std_region_data(KeyType const key)
{
  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (m_cached_entry->first == key))
  {
    return PtrHandle<T>(m_cached_entry->second);
  }

  // If the attempt to reuse the last entry failed, search for the key in the
  // map
  raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.cend())
  {
    m_cached_entry = it;
    return PtrHandle<T>(it->second);
  }

  m_cached_entry = m_std_region_data.cend();
  return PtrHandle<T>();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
const PtrHandle<T const> DataMapImplNoThreads<KeyType, T, Hasher>::std_region_data(
    KeyType const key) const
{
  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (m_cached_entry->first == key))
  {
    return PtrHandle<T const>(m_cached_entry->second);
  }

  const_raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T const>(it->second);
  }

  return PtrHandle<T const>();
}

// ----------------------------------------------------------------------------
// Implementation using threads
// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher = typename KeyType::hash_type>
class DataMapImplWithThreads
{
  private:
  using const_iterator_impl = typename std::unordered_map<KeyType, T *, Hasher>::const_iterator;
  using iterator_impl       = typename std::unordered_map<KeyType, T *, Hasher>::iterator;

  public:
  /// TYPEDEFS

  using key_type = KeyType;

  using const_iterator = DataMapIterator<T, DataMapConstItTraits<T, KeyType, const_iterator_impl>>;

  using iterator = DataMapIterator<T, DataMapItTraits<T, KeyType, iterator_impl>>;

  /// METHODS

  /// Default constructor
  DataMapImplWithThreads();

  /// Default destructor
  ~DataMapImplWithThreads();

  /// Return size of the underlying map
  Uint size() const;

  /// Remove one entry
  // void remove(KeyType const key);

  /// Remove all internal data
  void clear();

  /// Create a new entry
  PtrHandle<T> create(KeyType const key);

  /// Create a new entry, const version
  const PtrHandle<T const> create(KeyType const key) const;

  /// Get one entry
  PtrHandle<T> std_region_data(KeyType const key);

  /// Get one entry, const version
  const PtrHandle<T const> std_region_data(KeyType const key) const;

  /// Return the constant 'begin' iterator over the data
  inline const_iterator cbegin() const
  {
    const_iterator it(m_std_region_data.cbegin());
    return it;
  }

  /// Return the constant 'end' iterator over the data
  inline const_iterator cend() const
  {
    const_iterator it(m_std_region_data.cend());
    return it;
  }

  inline iterator begin()
  {
    iterator it(m_std_region_data.begin());
    return it;
  }

  inline iterator end()
  {
    iterator it(m_std_region_data.end());
    return it;
  }

  private:
  /// TYPES
  /// To use a custom type as key in unordered map, two things have to be
  /// defined: 1) An equality operator:
  ///    const bool operator==(const Key &other) const
  ///    {
  ///      ...
  ///    }
  /// 2) A hash function (see DefaultStdRegMapKeyHasher at the top of this
  /// file)

  using raw_iterator       = typename std::unordered_map<KeyType, T *, Hasher>::iterator;
  using const_raw_iterator = typename std::unordered_map<KeyType, T *, Hasher>::const_iterator;

  /// Map [std. region id => data for this type of std region]
  std::unordered_map<KeyType, T *, Hasher> m_std_region_data;
  const_raw_iterator m_cached_entry;

  /// Mutex for locking when map is modified
  mutable std::mutex m_mutex;
};

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
DataMapImplWithThreads<KeyType, T, Hasher>::DataMapImplWithThreads()
{
  m_cached_entry = m_std_region_data.end();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
DataMapImplWithThreads<KeyType, T, Hasher>::~DataMapImplWithThreads()
{
  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete it->second;
    it->second = nullptr;
  }
  m_std_region_data.clear();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
Uint DataMapImplWithThreads<KeyType, T, Hasher>::size() const
{
  std::lock_guard<std::mutex> lock(m_mutex);
  return m_std_region_data.size();
}

// ----------------------------------------------------------------------------

/*
template <typename KeyType, typename T, typename Hasher>
void StdRegionDataMapImplWithThreads<KeyType, T, Hasher>::remove(KeyType const
key)
{
  std::lock_guard<std::mutex> lock(m_mutex);

  raw_iterator it = m_std_region_data.find(key);
  if (it != m_std_region_data.end())
  {
    delete it->second;
    m_std_region_data.erase(it);
  }

  if (!m_std_region_data.empty())
  {
    m_cached_entry = m_std_region_data.cbegin();
  }
  else
  {
    m_cached_entry = m_std_region_data.cend();
  }
}
*/

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
void DataMapImplWithThreads<KeyType, T, Hasher>::clear()
{
  std::lock_guard<std::mutex> lock(m_mutex);

  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete it->second;
  }
  m_std_region_data.clear();
  m_cached_entry = m_std_region_data.cend();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
PtrHandle<T> DataMapImplWithThreads<KeyType, T, Hasher>::create(KeyType const key)
{
  std::lock_guard<std::mutex> lock(m_mutex);

  raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T>(it->second);
  }

  T *new_data = new T();
  m_std_region_data.insert(std::pair<KeyType, T *>(key, new_data));
  return PtrHandle<T>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
const PtrHandle<T const> DataMapImplWithThreads<KeyType, T, Hasher>::create(KeyType const key) const
{
  std::lock_guard<std::mutex> lock(m_mutex);

  const_raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T const>(it->second);
  }

  T *new_data = new T();
  m_std_region_data.insert(std::pair<KeyType, T *>(key, new_data));
  return PtrHandle<T const>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
PtrHandle<T> DataMapImplWithThreads<KeyType, T, Hasher>::std_region_data(KeyType const key)
{
  std::lock_guard<std::mutex> lock(m_mutex);

  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (m_cached_entry->first == key))
  {
    return PtrHandle<T>(m_cached_entry->second);
  }

  // If the attempt to reuse the last entry failed, search for the key in the
  // map
  raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.cend())
  {
    m_cached_entry = it;
    return PtrHandle<T>(it->second);
  }

  m_cached_entry = m_std_region_data.cend();
  return PtrHandle<T>();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T, typename Hasher>
const PtrHandle<T const> DataMapImplWithThreads<KeyType, T, Hasher>::std_region_data(
    KeyType const key) const
{
  std::lock_guard<std::mutex> lock(m_mutex);

  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (m_cached_entry->first == key))
  {
    return PtrHandle<T const>(m_cached_entry->second);
  }

  const_raw_iterator it = m_std_region_data.find(key);

  if (it != m_std_region_data.end())
  {
    return PtrHandle<T const>(it->second);
  }

  return PtrHandle<T const>();
}

// ----------------------------------------------------------------------------
// Implementation using std::vector, no threads
// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
class VecDataMapImplNoThreads
{
  private:
  using const_iterator_impl = typename std::vector<std::tuple<KeyType, T *>>::const_iterator;
  using iterator_impl       = typename std::vector<std::tuple<KeyType, T *>>::iterator;

  public:
  /// TYPEDEFS

  using key_type = KeyType;

  using const_iterator = DataMapIterator<T, DataMapConstItTraits<T, KeyType, const_iterator_impl>>;

  using iterator = DataMapIterator<T, DataMapItTraits<T, KeyType, iterator_impl>>;

  /// METHODS

  /// Default constructor
  VecDataMapImplNoThreads();

  /// Default destructor
  ~VecDataMapImplNoThreads();

  /// Return size of the underlying map
  Uint size() const;

  /// Remove one entry
  // void remove(KeyType const key);

  /// Remove all internal data
  void clear();

  /// Create a new entry
  PtrHandle<T> create(KeyType const key);

  /// Create a new entry, const version
  const PtrHandle<T const> create(KeyType const key) const;

  /// Get one entry
  PtrHandle<T> std_region_data(KeyType const key);

  /// Get one entry, const version
  const PtrHandle<T const> std_region_data(KeyType const key) const;

  /// Return the constant 'begin' iterator over the data
  inline const_iterator cbegin() const
  {
    const_iterator it(m_std_region_data.cbegin());
    return it;
  }

  /// Return the constant 'end' iterator over the data
  inline const_iterator cend() const
  {
    const_iterator it(m_std_region_data.cend());
    return it;
  }

  inline iterator begin()
  {
    iterator it(m_std_region_data.begin());
    return it;
  }

  inline iterator end()
  {
    iterator it(m_std_region_data.end());
    return it;
  }

  private:
  /// TYPES
  using raw_iterator       = typename std::vector<std::tuple<KeyType, T *>>::iterator;
  using const_raw_iterator = typename std::vector<std::tuple<KeyType, T *>>::const_iterator;

  std::vector<std::tuple<KeyType, T *>> m_std_region_data;
  const_raw_iterator m_cached_entry;
};

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
VecDataMapImplNoThreads<KeyType, T>::VecDataMapImplNoThreads()
{
  m_cached_entry = m_std_region_data.cend();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
VecDataMapImplNoThreads<KeyType, T>::~VecDataMapImplNoThreads()
{
  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete std::get<1>(*it);
    std::get<1>(*it) = nullptr;
  }
  m_std_region_data.resize(0);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
Uint VecDataMapImplNoThreads<KeyType, T>::size() const
{
  return m_std_region_data.size();
}

// ----------------------------------------------------------------------------

/*
template <typename KeyType, typename T>
void StdRegionVecMapImplNoThreads<KeyType, T>::remove(KeyType const key)
{
  Uint remove_pos = m_std_region_data.size();

  for (Uint i = 0; i < m_std_region_data.size(); ++i)
  {
    if (std::get<0>(m_std_region_data[i]) == key)
    {
      remove_pos = i;
      break;
    }
  }

  if (remove_pos < m_std_region_data.size())
  {
    delete std::get<1>(m_std_region_data[remove_pos]);

    if ((remove_pos + 1) < m_std_region_data.size())
    {
      std::swap(m_std_region_data[remove_pos], m_std_region_data.back());
    }
    m_std_region_data.pop_back();
  }

  if (!m_std_region_data.empty())
  {
    m_cached_entry = m_std_region_data.cbegin();
  }
  else
  {
    m_cached_entry = m_std_region_data.cend();
  }
}
*/

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
void VecDataMapImplNoThreads<KeyType, T>::clear()
{
  for (raw_iterator it = m_std_region_data.begin(); it != m_std_region_data.end(); ++it)
  {
    delete std::get<1>(*it);
    std::get<1>(*it) = nullptr;
  }
  m_std_region_data.resize(0);
  m_cached_entry = m_std_region_data.cend();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
PtrHandle<T> VecDataMapImplNoThreads<KeyType, T>::create(KeyType const key)
{
  for (Uint i = 0; i < m_std_region_data.size(); ++i)
  {
    if (std::get<0>(m_std_region_data[i]) == key)
    {
      return PtrHandle<T>(std::get<1>(m_std_region_data[i]));
    }
  }

  T *new_data = new T();
  m_std_region_data.push_back(std::tuple<KeyType, T *>(key, new_data));

  // Reset cached entry iterator, because pushing back to vector might
  // reallocate it, thus invalidating iterators
  m_cached_entry = m_std_region_data.cend();
  return PtrHandle<T>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
const PtrHandle<T const> VecDataMapImplNoThreads<KeyType, T>::create(KeyType const key) const
{
  for (Uint i = 0; i < m_std_region_data.size(); ++i)
  {
    if (std::get<0>(m_std_region_data[i]) == key)
    {
      return PtrHandle<const T>(std::get<1>(m_std_region_data[i]));
    }
  }

  T *new_data = new T();
  m_std_region_data.push_back(std::tuple<KeyType, T *>(key, new_data));

  // Reset cached entry iterator, because pushing back to vector might
  // reallocate it, thus invalidating iterators
  m_cached_entry = m_std_region_data.cend();
  return PtrHandle<T const>(new_data);
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
PtrHandle<T> VecDataMapImplNoThreads<KeyType, T>::std_region_data(KeyType const key)
{
  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (std::get<0>(*m_cached_entry) == key))
  {
    return PtrHandle<T>(std::get<1>(*m_cached_entry));
  }

  // If the attempt to reuse the last entry failed, search for the key in the
  // map
  for (Uint i = 0; i < m_std_region_data.size(); ++i)
  {
    if (std::get<0>(m_std_region_data[i]) == key)
    {
      m_cached_entry = m_std_region_data.begin() + i;
      return PtrHandle<T>(std::get<1>(m_std_region_data[i]));
    }
  }

  m_cached_entry = m_std_region_data.cend();
  return PtrHandle<T>();
}

// ----------------------------------------------------------------------------

template <typename KeyType, typename T>
const PtrHandle<T const> VecDataMapImplNoThreads<KeyType, T>::std_region_data(
    KeyType const key) const
{
  // Try to reuse the last entry without searching the map
  if ((m_cached_entry != m_std_region_data.cend()) && (std::get<0>(*m_cached_entry) == key))
  {
    return PtrHandle<const T>(std::get<1>(*m_cached_entry));
  }

  // If the attempt to reuse the last entry failed, search for the key in the
  // map
  for (Uint i = 0; i < m_std_region_data.size(); ++i)
  {
    if (std::get<0>(m_std_region_data[i]) == key)
    {
      // m_cached_entry = m_std_region_data.cbegin() + i;
      return PtrHandle<const T>(std::get<1>(m_std_region_data[i]));
    }
  }

  // m_cached_entry = m_std_region_data.cend();
  return PtrHandle<const T>();
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace common

} // namespace pdekit

#endif
