#ifndef PDEKIT_Common_FW_Instance_Container_hpp
#define PDEKIT_Common_FW_Instance_Container_hpp

#include <map>
#include <mutex>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

namespace detail
{
template <typename T, typename FwPolicy>
class FWInstanceContainer
{
  public:
  using key_type       = typename FwPolicy::key_type;
  using value_type     = T;
  using catalog_type   = std::map<typename FwPolicy::key_type, T>;
  using iterator       = typename catalog_type::iterator;
  using const_iterator = typename catalog_type::const_iterator;

  /// Default constructor
  FWInstanceContainer();

  /// Copy constructor - delete, because std::mutex has deleted
  /// copy constructor anyway
  FWInstanceContainer(const FWInstanceContainer &other) = delete;

  /// Delete assignment operator
  FWInstanceContainer &operator=(const FWInstanceContainer &rhs) = delete;

  /// Default destructor
  ~FWInstanceContainer();

  /// Set iterator to point to an instance based on key
  const_iterator find_instance(const key_type &key) const;

  /// Create a new instance
  const_iterator create_instance(const key_type &key);

  /// Iterator to the beginning
  // iterator begin() noexcept;

  /// Const iterator to the beginning
  const_iterator cbegin() const noexcept;

  /// Iterator to the end
  // iterator end() noexcept;

  /// Const iterator to the end
  const_iterator cend() const noexcept;

  private:
  /// Container for flyweight instances
  catalog_type m_catalog;

  /// Mutex for locking when new instance is created
  std::mutex m_mutex;
};

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
FWInstanceContainer<T, FwPolicy>::FWInstanceContainer()
{
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
FWInstanceContainer<T, FwPolicy>::~FWInstanceContainer()
{
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::const_iterator FWInstanceContainer<
    T, FwPolicy>::find_instance(const key_type &key) const
{
  const_iterator it = m_catalog.find(key);
  return it;
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::const_iterator FWInstanceContainer<
    T, FwPolicy>::create_instance(const key_type &key)
{
  std::lock_guard<std::mutex> lock(m_mutex);

  iterator it = m_catalog.find(key);

  if (it == m_catalog.end())
  {
    m_catalog.insert(std::pair<key_type, T>(key, T()));
    it = m_catalog.find(key);
    T::construct(key, it->second);
  }

  // We're returning an iterator, which should be implicitly convertible to
  // const_iterator (in STL containers at least)
  return it;
}

// ----------------------------------------------------------------------------

/*
template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::iterator
FWInstanceContainer<T, FwPolicy>::begin() noexcept
{
  return m_catalog.begin();
}
*/

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::const_iterator FWInstanceContainer<T, FwPolicy>::cbegin()
    const noexcept
{
  return m_catalog.cbegin();
}

// ----------------------------------------------------------------------------

/*
template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::iterator FWInstanceContainer<T,
FwPolicy>::end() noexcept
{
  return m_catalog.end();
}
*/

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
typename FWInstanceContainer<T, FwPolicy>::const_iterator FWInstanceContainer<T, FwPolicy>::cend()
    const noexcept
{
  return m_catalog.cend();
}

// ----------------------------------------------------------------------------

} // namespace detail

} // namespace common

} // namespace pdekit

#endif
