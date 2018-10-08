#ifndef PDEKIT_Common_Data_Map_Iterator_hpp
#define PDEKIT_Common_Data_Map_Iterator_hpp

#include "common/Constants.hpp"
#include "common/PtrHandle.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
class DataMapIterator
{
  /// TYPEDEFS
  public:
  using ptr_handle_type = PtrHandle<typename Traits::value_type>;
  using key_type        = typename Traits::key_type;
  using reference_type  = typename Traits::reference_type;

  private:
  using wrapped_iterator_type = typename Traits::wrapped_iterator_type;

  public:
  /// Default constructor
  DataMapIterator();

  /// Construct from wrapped iterator
  DataMapIterator(wrapped_iterator_type it);

  /// Increment the iterator
  inline DataMapIterator &operator++();

  /// Decrease the operator
  inline DataMapIterator &operator--();

  /// Get the standard region tag of the entry that
  /// the iterator points to
  inline key_type key_value() const;

  /// Data of the entry to which the iterator points to
  inline ptr_handle_type data_ptr();

  /// Compare two iterators and return true if they are equal
  inline bool operator==(const DataMapIterator &rhs_it) const
  {
    return m_iterator == rhs_it.m_iterator;
  }

  /// Compare two iterators and return true if they are different
  inline bool operator!=(const DataMapIterator &rhs_it) const
  {
    return m_iterator != rhs_it.m_iterator;
  }

  /// Get the value the iterator points to
  reference_type operator*();

  /// Get the address (pointer) the iterator points to
  ptr_handle_type operator->();

  private:
  /// MEMBER VARIABLES

  /// The wrapped iterator - this is the actual iterator over the std region
  /// data map that we increment/decrement
  wrapped_iterator_type m_iterator;
};

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
DataMapIterator<T, Traits>::DataMapIterator()
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
DataMapIterator<T, Traits>::DataMapIterator(wrapped_iterator_type it) : m_iterator(it)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
DataMapIterator<T, Traits> &DataMapIterator<T, Traits>::operator++()
{
  m_iterator++;
  return *this;
}
// ----------------------------------------------------------------------------

template <typename T, typename Traits>
DataMapIterator<T, Traits> &DataMapIterator<T, Traits>::operator--()
{
  m_iterator--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
typename DataMapIterator<T, Traits>::key_type DataMapIterator<T, Traits>::key_value() const
{
  // return m_iterator->first;
  return Traits::get_key_from_iterator(m_iterator);
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
typename DataMapIterator<T, Traits>::ptr_handle_type DataMapIterator<T, Traits>::data_ptr()
{
  // ptr_handle_type ptr(m_iterator->second);
  ptr_handle_type ptr(Traits::get_value_ptr_from_iterator(m_iterator));
  return ptr;
}

// ----------------------------------------------------------------------------

/*
template <typename T, typename Traits>
inline bool operator==(const DataMapIterator<T, Traits> &it_lhs,
                             const DataMapIterator<T, Traits> &it_rhs)
{
  return it_lhs.m_iterator == it_rhs.m_iterator;
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
inline bool operator!=(const DataMapIterator<T, Traits> &it_lhs,
                             const DataMapIterator<T, Traits> &it_rhs)
{
  return it_lhs.m_iterator != it_rhs.m_iterator;
}
*/

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
typename DataMapIterator<T, Traits>::reference_type DataMapIterator<T, Traits>::operator*()
{
  // ptr_handle_type ptr_handle(m_iterator->second);
  ptr_handle_type ptr_handle(Traits::get_value_ptr_from_iterator(m_iterator));
  return *ptr_handle;
}

// ----------------------------------------------------------------------------

template <typename T, typename Traits>
typename DataMapIterator<T, Traits>::ptr_handle_type DataMapIterator<T, Traits>::operator->()
{
  // ptr_handle_type ptr_handle(m_iterator->second);
  ptr_handle_type ptr_handle(Traits::get_value_ptr_from_iterator(m_iterator));
  return ptr_handle;
}

// ----------------------------------------------------------------------------
// Traits for const and non-const iterator over standard region data map
// ----------------------------------------------------------------------------

template <typename T, typename KeyType, typename IteratorImpl>
struct DataMapConstItTraits
{
  using value_type            = T const;
  using key_type              = KeyType;
  using reference_type        = T const &;
  using pointer_type          = T const *;
  using wrapped_iterator_type = IteratorImpl;

  inline static key_type get_key_from_iterator(const IteratorImpl &it)
  {
    return std::get<0>(*it);
  }

  inline static pointer_type get_value_ptr_from_iterator(const IteratorImpl &it)
  {
    return std::get<1>(*it);
  }

  /*
  inline static key_type
  get_key_from_iterator(const typename std::vector<std::tuple<KeyType, T
  *>>::const_iterator &it)
  {
    return std::get<0>(*it);
  }

  inline static pointer_type get_value_ptr_from_iterator(
      const typename std::vector<std::tuple<KeyType, T *>>::const_iterator
  &it)
  {
    return std::get<1>(*it);
  }
  */
};

template <typename T, typename KeyType, typename IteratorImpl>
struct DataMapItTraits
{
  using value_type            = T;
  using key_type              = KeyType;
  using reference_type        = T &;
  using pointer_type          = T *;
  using wrapped_iterator_type = IteratorImpl;

  inline static key_type get_key_from_iterator(const IteratorImpl &it)
  {
    return std::get<0>(*it);
  }

  inline static pointer_type get_value_ptr_from_iterator(const IteratorImpl &it)
  {
    return std::get<1>(*it);
  }

  /*
  inline static key_type
  get_key_from_iterator(const typename std::vector<std::tuple<KeyType, T
  *>>::const_iterator &it)
  {
    return std::get<0>(*it);
  }

  inline static pointer_type get_value_ptr_from_iterator(
      const typename std::vector<std::tuple<KeyType, T *>>::const_iterator
  &it)
  {
    return std::get<1>(*it);
  }
  */
};

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
