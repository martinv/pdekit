#ifndef PDEKIT_Common_Cell_Iterator_Range_hpp
#define PDEKIT_Common_Cell_Iterator_Range_hpp

#include <cstddef>

namespace pdekit
{

namespace common
{

template <typename IterType>
class IteratorRange
{
  public:
  /// TYPEDEFS

  using iterator       = IterType;
  using const_iterator = IterType;

  /// Default constructor: empty range
  IteratorRange();

  /// Constructor: take the first and last iterator
  /// as inputs
  IteratorRange(const IterType &first, const IterType &last);

  /// Copy constructor
  IteratorRange(const IteratorRange &other_range);

  /// Assignement operator
  IteratorRange &operator=(const IteratorRange &other_range);

  /// Destructor
  ~IteratorRange();

  /// Return the first iterator - beginning of the cell range
  IterType begin() const;

  /// Return the last iterator - one after the end of the cell range
  IterType end() const;

  /// Get the distance (length of the range)
  std::ptrdiff_t distance() const;

  private:
  /// Pointer/iterator to the first and last element in the range
  IterType m_begin;
  IterType m_end;
  std::ptrdiff_t m_distance;
};

// ----------------------------------------------------------------------------

template <typename IterType>
IteratorRange<IterType>::IteratorRange()
{
}

// ----------------------------------------------------------------------------

template <typename IterType>
IteratorRange<IterType>::IteratorRange(const IterType &first, const IterType &last)
    : m_begin(first), m_end(last), m_distance(IterType::distance(first, last))
{
}

// ----------------------------------------------------------------------------

template <typename IterType>
IteratorRange<IterType>::IteratorRange(const IteratorRange &other_range)
    : m_begin(other_range.m_begin), m_end(other_range.m_end), m_distance(other_range.m_distance)
{
}

// ----------------------------------------------------------------------------

template <typename IterType>
IteratorRange<IterType> &IteratorRange<IterType>::operator=(
    const IteratorRange<IterType> &other_range)
{
  m_begin    = other_range.m_begin;
  m_end      = other_range.m_end;
  m_distance = other_range.m_distance;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename IterType>
IteratorRange<IterType>::~IteratorRange()
{
}

// ----------------------------------------------------------------------------

template <typename IterType>
inline IterType IteratorRange<IterType>::begin() const
{
  return m_begin;
}

// ----------------------------------------------------------------------------

template <typename IterType>
inline IterType IteratorRange<IterType>::end() const
{
  return m_end;
}

// ----------------------------------------------------------------------------

template <typename IterType>
inline std::ptrdiff_t IteratorRange<IterType>::distance() const
{
  return m_distance;
}

// ----------------------------------------------------------------------------

// Free function to create an iterator range
template <typename IterType>
IteratorRange<IterType> make_iter_range(const IterType &begin, const IterType &end)
{
  return IteratorRange<IterType>(begin, end);
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
