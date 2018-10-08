#ifndef PDEKIT_Common_Range1D_hpp
#define PDEKIT_Common_Range1D_hpp

#include <initializer_list>
#include <ostream>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

template <typename IntType>
class Range1D
{
  public:
  /// Default constructor
  Range1D();

  /// Constructor from lower and upper bound
  explicit Range1D(const IntType lower, const IntType upper);

  /// Construct from initializer list
  Range1D(std::initializer_list<IntType> list);

  /// Destructor
  ~Range1D();

  /// Get lower bound
  IntType lbound() const;

  /// Get upper bound
  IntType ubound() const;

  /// \brief Return the size of the range (ubound() - lbound() + 1)
  IntType size() const;

  /// \brief Return true if the range is empty, i.e. the lower bound
  /// is larger than upper bound
  bool empty() const;

  /// \brief Return true if the index is in range
  bool in_range(IntType i) const;

  /// \brief Increment the range by a constant
  Range1D &operator+=(IntType incr);

  /// \brief Deincrement the range by a constant
  Range1D &operator-=(IntType incr);

  private:
  IntType m_lower;
  IntType m_upper;
};

// ----------------------------------------------------------------------------

template <typename IntType>
Range1D<IntType>::Range1D() : m_lower(IntType()), m_upper(IntType())
{
}

// ----------------------------------------------------------------------------

template <typename IntType>
Range1D<IntType>::Range1D(const IntType lower, const IntType upper) : m_lower(lower), m_upper(upper)
{
}

// ----------------------------------------------------------------------------

template <typename IntType>
Range1D<IntType>::Range1D(std::initializer_list<IntType> list)
{
  typename std::initializer_list<IntType>::const_iterator it = list.begin();

  if (list.size() > 0)
  {
    m_lower = *it;
  }
  it++;
  if (it != list.end())
  {
    m_upper = *it;
  }
}

// ----------------------------------------------------------------------------

template <typename IntType>
Range1D<IntType>::~Range1D()
{
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline IntType Range1D<IntType>::lbound() const
{
  return m_lower;
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline IntType Range1D<IntType>::ubound() const
{
  return m_upper;
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline IntType Range1D<IntType>::size() const
{
  return m_upper - m_lower + 1;
}

// ----------------------------------------------------------------------------

template <typename IntType>
bool Range1D<IntType>::empty() const
{
  return m_upper < m_lower;
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline bool Range1D<IntType>::in_range(IntType i) const
{
  return (m_lower <= i) && (i <= m_upper);
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline Range1D<IntType> &Range1D<IntType>::operator+=(IntType incr)
{
  m_lower += incr;
  m_upper += incr;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename IntType>
inline Range1D<IntType> &Range1D<IntType>::operator-=(IntType incr)
{
  m_lower -= incr;
  m_upper -= incr;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename IntType>
std::ostream &operator<<(std::ostream &os, const Range1D<IntType> &range)
{
  os << "[" << range.lbound() << "," << range.ubound() << "]";
  return os;
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif // PDEKIT_Common_Range1D_hpp
