#ifndef PDEKIT_Common_Array_View_hpp
#define PDEKIT_Common_Array_View_hpp

#include <ostream>
#include <vector>

#include "common/Constants.hpp"

namespace pdekit
{

namespace common
{

/// ----------------------------------------------------------------------------
/// ArrayView is used as a lightweight proxy to represent one row of another
/// class which holds a table of values. For example, it is used to obtain
/// one data block from BlockArray
/// The class ArrayView itself does not hold the row data, only knows where to
/// look for them. It is not responsible for destruction of the underlying data
/// to which it provides access.
/// We suppose that the 'Table' class stores each row as a continuous array
/// of values
/// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Case when the table slice is not constant, i.e.
// the values it refers to CAN be modified
// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
class ArrayView
{
  public:
  /// TYPEDEFS
  using value_type      = T;
  using iterator        = typename std::conditional<std::is_const<T>::value, const T *, T *>::type;
  using const_iterator  = const T *;
  using reference       = typename std::conditional<std::is_const<T>::value, const T &, T &>::type;
  using const_reference = const T &;
  using size_type       = SizeType;
  using difference_type = std::ptrdiff_t;

  /// Default constructor
  ArrayView();

  /// Construct one row given pointer to first element and size
  /// The parameter @param first is not const, because the member variable
  ///                      m_first has to remain non-const so that we can
  ///                      modify the underlying values that m_first points to
  ArrayView(value_type *first, const size_type size);

  /// Copy constructor
  ArrayView(const ArrayView &other_slice);

  /// Assignement operator
  ArrayView &operator=(const ArrayView &other_slice);

  /// Destructor
  ~ArrayView();

  /// access to one element, read-only
  inline const_reference operator[](const size_type i) const;

  /// access to one element, writable
  inline reference operator[](const size_type i);

  /// Return the length of this table
  const size_type size() const;

  /// front() and back()
  reference front();
  const_reference front() const;

  reference back();
  const_reference back() const;

  private:
  /// Pointer to the first value
  value_type *m_first;

  /// Length of this row
  size_type m_size;
};

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
ArrayView<T, Dim, SizeType>::ArrayView() : m_first(nullptr), m_size(0)
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
ArrayView<T, Dim, SizeType>::ArrayView(value_type *first, const size_type size)
    : m_first(first), m_size(size)
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
ArrayView<T, Dim, SizeType>::ArrayView(const ArrayView &other_slice)
    : m_first(other_slice.m_first), m_size(other_slice.m_size)
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
ArrayView<T, Dim, SizeType> &ArrayView<T, Dim, SizeType>::operator=(const ArrayView &other_slice)
{
  m_first = other_slice.m_first;
  m_size  = other_slice.m_size;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
ArrayView<T, Dim, SizeType>::~ArrayView()
{
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::const_reference ArrayView<T, Dim, SizeType>::operator[](
    const size_type i) const
{
  return m_first[i];
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::reference ArrayView<T, Dim, SizeType>::operator[](
    const size_type i)
{
  return m_first[i];
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
const typename ArrayView<T, Dim, SizeType>::size_type ArrayView<T, Dim, SizeType>::size() const
{
  return m_size;
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::reference ArrayView<T, Dim, SizeType>::front()
{
  return *m_first;
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::const_reference ArrayView<T, Dim, SizeType>::front() const
{
  return *m_first;
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::reference ArrayView<T, Dim, SizeType>::back()
{
  return m_first[m_size - 1];
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
typename ArrayView<T, Dim, SizeType>::const_reference ArrayView<T, Dim, SizeType>::back() const
{
  return m_first[m_size - 1];
}

// ----------------------------------------------------------------------------

template <typename T, typename SizeType>
ArrayView<const T, _1D, SizeType> make_view(const std::vector<T> &vec)
{
  return ArrayView<const T, _1D, SizeType>(vec.data(), static_cast<SizeType>(vec.size()));
}

// ----------------------------------------------------------------------------

template <typename T, Uint Dim, typename SizeType>
std::ostream &operator<<(std::ostream &os, const ArrayView<T, Dim, SizeType> &slice)
{

  for (SizeType i = 0; i < (slice.size() - 1); ++i)
  {
    os << slice[i] << " ";
  }

  os << slice.back();
  return os;
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_Array_View_hpp
