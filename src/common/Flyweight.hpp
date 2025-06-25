#ifndef PDEKIT_Common_Flyweight_hpp
#define PDEKIT_Common_Flyweight_hpp

#include "common/FWInstanceContainer.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

struct DefaultFlyweightPolicy
{
  using key_type = Uint;
};

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy = DefaultFlyweightPolicy>
class Flyweight
{
  public:
  using catalog_type = detail::FWInstanceContainer<T, FwPolicy>;
  using key_type     = typename catalog_type::key_type;
  using value_type   = typename catalog_type::value_type;

  /// Default constructor
  Flyweight();

  /// Construct directly based on key
  Flyweight(const key_type &key);

  /// Copy constructor
  Flyweight(const Flyweight &other_fw);

  /// Assignment operator
  Flyweight &operator=(const Flyweight &other_fw);

  /// Destructor
  ~Flyweight();

  /// Get the key to which the current state of Flyweight
  /// corresponds
  const key_type &key() const;

  /// Get the underlying object
  const T &get() const;

  /// Return the number of instances so far created
  Uint nb_instances() const;

  /// Change the type of this flyweight
  template <typename... Params>
  void change_type(Params &&... params);

  private:
  /// Here we store instances
  static catalog_type Catalog;

  /// Iterator pointing to one instance of flyweight
  typename catalog_type::const_iterator m_instance;
};

// ----------------------------------------------------------------------------
// Initialize the static catalog:

template <typename T, typename FwPolicy>
typename Flyweight<T, FwPolicy>::catalog_type Flyweight<T, FwPolicy>::Catalog;

// Can't use assignment in initialization of static variable, because Catalog
// contains an instance of std::mutex, whose copy-constructor is not available
// (deleted)! template <typename T, typename FwPolicy> typename Flyweight<T,
// FwPolicy>::catalog_type
//     Flyweight<T, FwPolicy>::Catalog = Flyweight<T, FwPolicy>::catalog_type();
//

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Flyweight<T, FwPolicy>::Flyweight() : m_instance(nullptr)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Flyweight<T, FwPolicy>::Flyweight(const key_type &key)
{
  // change_type(key);
  if (key == T::undefined)
  {
    m_instance = typename catalog_type::iterator(nullptr);
    return;
  }

  m_instance = Catalog.find_instance(key);
  if (m_instance == Catalog.cend())
  {
    m_instance = Catalog.create_instance(key);
  }
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Flyweight<T, FwPolicy>::Flyweight(const Flyweight &other_fw) : m_instance(other_fw.m_instance)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Flyweight<T, FwPolicy> &Flyweight<T, FwPolicy>::operator=(const Flyweight &other_fw)
{
  m_instance = other_fw.m_instance;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Flyweight<T, FwPolicy>::~Flyweight<T, FwPolicy>()
{
  m_instance = typename catalog_type::iterator(nullptr);
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
inline const typename Flyweight<T, FwPolicy>::key_type &Flyweight<T, FwPolicy>::key() const
{
  return m_instance->first;
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
inline const T &Flyweight<T, FwPolicy>::get() const
{
  return m_instance->second;
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
Uint Flyweight<T, FwPolicy>::nb_instances() const
{
  return Catalog.size();
}

// ----------------------------------------------------------------------------

template <typename T, typename FwPolicy>
template <typename... Params>
void Flyweight<T, FwPolicy>::change_type(Params &&... params)
{
  key_type const key(std::forward<Params>(params)...);

  m_instance = Catalog.find_instance(key);
  if (m_instance == Catalog.cend())
  {
    m_instance = Catalog.create_instance(key);
  }
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
