#ifndef PDEKIT_Mesh_Key_Cache_hpp
#define PDEKIT_Mesh_Key_Cache_hpp

#include <vector>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace mesh
{

// ----------------------------------------------------------------------------

template <typename KeyType>
class KeyCache
{
  public:
  KeyCache() = default;

  /// Disable copy constructor
  KeyCache(const KeyCache &rhs) = delete;

  /// Default destructor
  ~KeyCache() = default;

  /// Disable assignment operator
  KeyCache &operator=(const KeyCache &rhs) = delete;

  /// Get the size of the cache
  Uint size() const;

  /// Clear the contents of the cache
  void clear();

  /// Reserve the size of cache
  void reserve(const Uint nb_entries);

  /// Add a new key
  void push_back(const KeyType key);

  /// Change the value of key
  void change_key(const Uint idx, const KeyType key);

  /// Return the value of key on position idx
  KeyType key(const Uint idx) const;

  private:
  Uint get_key_storage_idx(const KeyType key) const;

  /// Vector of keys
  std::vector<KeyType> m_keys;

  /// Vector that maps index to discrete element key
  std::vector<Uint> m_index_to_key_map;
};

// ----------------------------------------------------------------------------

template <typename KeyType>
Uint KeyCache<KeyType>::size() const
{
  return m_index_to_key_map.size();
}

// ----------------------------------------------------------------------------

template <typename KeyType>
void KeyCache<KeyType>::clear()
{
  m_keys.clear();
  m_index_to_key_map.clear();
}

// ----------------------------------------------------------------------------

template <typename KeyType>
void KeyCache<KeyType>::reserve(const Uint nb_entries)
{
  m_keys.clear();
  m_index_to_key_map.reserve(nb_entries);
}

// ----------------------------------------------------------------------------

template <typename KeyType>
void KeyCache<KeyType>::push_back(const KeyType key)
{
  Uint key_idx = get_key_storage_idx(key);

  if (key_idx == m_keys.size())
  {
    m_keys.push_back(key);
    key_idx = m_keys.size() - 1;
  }

  m_index_to_key_map.push_back(key_idx);
}

// ----------------------------------------------------------------------------

template <typename KeyType>
void KeyCache<KeyType>::change_key(const Uint idx, const KeyType key)
{
  Uint key_idx = get_key_storage_idx(key);
  if (key_idx == m_keys.size())
  {
    m_keys.push_back(key);
    key_idx = m_keys.size() - 1;
  }
  m_index_to_key_map[idx] = key_idx;
}

// ----------------------------------------------------------------------------

template <typename KeyType>
inline KeyType KeyCache<KeyType>::key(const Uint idx) const
{
  return m_keys[m_index_to_key_map[idx]];
}

// ----------------------------------------------------------------------------

template <typename KeyType>
Uint KeyCache<KeyType>::get_key_storage_idx(const KeyType key) const
{
  for (Uint i = 0; i < m_keys.size(); ++i)
  {
    if (m_keys[i] == key)
    {
      return i;
    }
  }

  return m_keys.size();
}

// ----------------------------------------------------------------------------

} // namespace mesh

} // namespace pdekit

#endif
