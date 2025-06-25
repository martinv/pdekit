#ifndef PDEKIT_Common_Tagged_Int_hpp
#define PDEKIT_Common_Tagged_Int_hpp

#include <iosfwd>
#include <type_traits>

namespace pdekit
{

namespace common
{

template <typename T, typename Tag>
class TaggedInt
{
  public:
  using value_type = T;

  /// Default constructor
  TaggedInt();

  /// Constructor from value
  /// @param val ... value to construct from
  explicit TaggedInt(const T val);

  /// Constructor from another index
  /// @param other ... index to construct from
  TaggedInt(const TaggedInt &other);

  /// Constructor from another index
  /// @param other ... index to construct from
  TaggedInt(TaggedInt &&other);

  /// Destructor
  ~TaggedInt();

  /// Assignment operator from value
  TaggedInt &operator=(const T val);

  /// Assignment operator from another value
  TaggedInt &operator=(const TaggedInt &rhs);

  /// Move assignment from another value
  TaggedInt &operator=(TaggedInt &&rhs);

  /// Increment index (prefix)
  TaggedInt &operator++();

  /// Increment index (postfix)
  TaggedInt operator++(int);

  /// Decrement index (prefix)
  TaggedInt &operator--();

  /// Decrement index (postfix)
  TaggedInt operator--(int);

  /// Get value
  /// @return ... stored value
  T id() const;

  private:
  static_assert(std::is_integral<T>::value,
                "The first template parameter in TaggedInt<T, Tag> must be of "
                "integer type.");
  T m_id;
};

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag>::TaggedInt() : m_id(T())
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag>::TaggedInt(const T val) : m_id(val)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag>::TaggedInt(const TaggedInt &other) : m_id(other.m_id)
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag>::TaggedInt(TaggedInt &&other) : m_id(other.m_id)
{
  other.m_id = 0;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag>::~TaggedInt()
{
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> &TaggedInt<T, Tag>::operator=(const T val)
{
  m_id = val;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> &TaggedInt<T, Tag>::operator=(const TaggedInt &rhs)
{
  m_id = rhs.m_id;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> &TaggedInt<T, Tag>::operator=(TaggedInt &&rhs)
{
  m_id     = rhs.m_id;
  rhs.m_id = 0;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> &TaggedInt<T, Tag>::operator++()
{
  m_id++;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> TaggedInt<T, Tag>::operator++(int)
{
  TaggedInt temp = *this;
  ++*this;
  return temp;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> &TaggedInt<T, Tag>::operator--()
{
  m_id--;
  return *this;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
TaggedInt<T, Tag> TaggedInt<T, Tag>::operator--(int)
{
  TaggedInt temp = *this;
  --*this;
  return temp;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
inline T TaggedInt<T, Tag>::id() const
{
  return m_id;
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
inline bool operator==(const TaggedInt<T, Tag> &lhs, const TaggedInt<T, Tag> &rhs)
{
  return lhs.id() == rhs.id();
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
inline bool operator!=(const TaggedInt<T, Tag> &lhs, const TaggedInt<T, Tag> &rhs)
{
  return lhs.id() != rhs.id();
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
inline bool operator<(const TaggedInt<T, Tag> &lhs, const TaggedInt<T, Tag> &rhs)
{
  return lhs.id() < rhs.id();
}

// ----------------------------------------------------------------------------

template <typename T, typename Tag>
std::ostream &operator<<(std::ostream &os, const TaggedInt<T, Tag> &idx)
{
  os << idx.id();
  return os;
}

// ----------------------------------------------------------------------------

} // namespace common

} // namespace pdekit

#endif
