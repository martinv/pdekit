#ifndef PDEKIT_Common_Small_Integer_Pack_hpp
#define PDEKIT_Common_Small_Integer_Pack_hpp

#include <iostream>

#include "common/Meta.hpp"

namespace pdekit
{

namespace common
{

namespace detail
{

// ----------------------------------------------------------------------------

// Forward declaration:
template <Uint NbTermsToSum, Uint... Sizes>
struct PartialSum;

template <Uint I, Uint... Sizes>
struct PartialSum<0, I, Sizes...>
{
  enum
  {
    value = 0
  };
};

template <Uint I, Uint... Sizes>
struct PartialSum<1, I, Sizes...>
{
  enum
  {
    value = I
  };
};

template <Uint NbTermsToSum, Uint I, Uint... Sizes>
struct PartialSum<NbTermsToSum, I, Sizes...>
{
  enum
  {
    value = I + PartialSum<NbTermsToSum - 1, Sizes...>::value
  };
};

// ----------------------------------------------------------------------------

template <Uint... Sizes>
struct Offset;

template <Uint... Sizes>
struct Offset<0, Sizes...>
{
  enum
  {
    value = 0u
  };
};

template <Uint Idx, Uint... Sizes>
struct Offset<Idx, Sizes...>
{
  enum
  {
    value = PartialSum<Idx, Sizes...>::value
  };
};

// ----------------------------------------------------------------------------

template <Uint... Sizes>
struct Mask;

template <Uint I, Uint... Sizes>
struct Mask<I, Sizes...>
{
  enum
  {
    value = PowerOfTwo<Offset<I + 1, Sizes...>::value>::value -
            PowerOfTwo<Offset<I, Sizes...>::value>::value
  };
};

template <Uint... Sizes>
struct Mask<0, Sizes...>
{
  enum
  {
    value = PowerOfTwo<Offset<1, Sizes...>::value>::value - 1
  };
};

// ----------------------------------------------------------------------------

template <Uint bits>
struct Store;
template <>
struct Store<8>
{
  typedef uint8_t type;
};
template <>
struct Store<16>
{
  typedef uint16_t type;
};

template <>
struct Store<18>
{
  typedef uint32_t type;
};

template <>
struct Store<24>
{
  typedef uint32_t type;
};

template <>
struct Store<28>
{
  typedef uint32_t type;
};

template <>
struct Store<30>
{
  typedef uint32_t type;
};

template <>
struct Store<32>
{
  typedef uint32_t type;
};
template <>
struct Store<64>
{
  typedef uint64_t type;
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------
// A class that holds several (small) integer values by occupying different
// bits of one integer which is if of type 'Store< >'. Note that the store
// is specified only for 8, 16, 32 and 64 bits. This means that when
// SmallIntegerPack is instantiated, the sum of the template parameter
// called 'Sizes' must be one of these values:
// Example:
// SmallIntegerPack<1,2,5,8> sip1; // This is fine, sum is 16
// SmallIntegerPack<1,5> sip2;     // This will result in compilation
//                                 // error: sum is 6
// SmallIngeterPack<3,5> sip3;     // OK
// ----------------------------------------------------------------------------

template <Uint... Sizes>
class SmallIntegerPack
{
  public:
  /// TYPEDEFS
  typedef typename detail::Store<Sum<Sizes...>::value>::type store_type;

  enum
  {
    nb_fields      = sizeof...(Sizes),
    nb_stored_bits = Sum<Sizes...>::value
  };

  /// METHODS

  /// Default constructor
  SmallIntegerPack();

  /// Construct from value (across all bits)
  constexpr SmallIntegerPack(store_type value);

  /// Copy constructor
  constexpr SmallIntegerPack(const SmallIntegerPack &other);

  /// Assignment operator
  SmallIntegerPack &operator=(const SmallIntegerPack &other);

  /// Destructor
  ~SmallIntegerPack();

  /// Return the tag value
  constexpr Uint store_value() const
  {
    return m_store;
  }

  /// Return the value stored in field on position 'pos'
  template <Uint pos>
  constexpr Uint field() const
  {
    return ((m_store & detail::Mask<pos, Sizes...>::value) >> detail::Offset<pos, Sizes...>::value);
  }

  /// Set the field value on position 'pos'
  template <Uint pos>
  inline void set_field(const Uint field_value)
  {
    m_store &= ~detail::Mask<pos, Sizes...>::value;
    const Uint new_mask = (field_value << detail::Offset<pos, Sizes...>::value);
    m_store |= new_mask;
  }

  /// Compare two element types, check if they are the same
  inline bool operator==(const SmallIntegerPack &other) const
  {
    return (m_store == other.m_store);
  }

  /// Compare two element types, check if they are different
  inline bool operator!=(const SmallIntegerPack &other) const
  {
    return (m_store != other.m_store);
  }

  /// Comparison operator for ordering
  inline bool operator<(const SmallIntegerPack &other) const
  {
    return m_store < other.m_store;
  }

  /// Print for debugging
  void print() const
  {
    std::cout << "Stored value = " << m_store << std::endl;
  }

  /// Print the configuration data: offests and masks. This is needed for
  /// debugging
  void print_config() const
  {
    std::cout << "*****************************************" << std::endl;
    std::cout << "Configuration data" << std::endl;
    std::cout << "Offset 0 = " << detail::Offset<0, Sizes...>::value << std::endl;
    std::cout << "Offset 1 = " << detail::Offset<1, Sizes...>::value << std::endl;
    std::cout << "Offset 2 = " << detail::Offset<2, Sizes...>::value << std::endl;
    std::cout << "Offset 3 = " << detail::Offset<3, Sizes...>::value << std::endl;

    std::cout << "Mask 0 = " << detail::Mask<0, Sizes...>::value << std::endl;
    std::cout << "Mask 1 = " << detail::Mask<1, Sizes...>::value << std::endl;
    std::cout << "Mask 2 = " << detail::Mask<2, Sizes...>::value << std::endl;
    std::cout << "Mask 3 = " << detail::Mask<3, Sizes...>::value << std::endl;
    std::cout << "*****************************************" << std::endl;
  }

  protected:
  store_type m_store;

  // std::tuple<decltype(Sizes)...> offsets;
  //
};

// ----------------------------------------------------------------------------

template <Uint... Sizes>
SmallIntegerPack<Sizes...>::SmallIntegerPack() : m_store(store_type())
{
}

// ----------------------------------------------------------------------------

template <Uint... Sizes>
constexpr SmallIntegerPack<Sizes...>::SmallIntegerPack(store_type value) : m_store(value)
{
}

// ----------------------------------------------------------------------------

template <Uint... Sizes>
constexpr SmallIntegerPack<Sizes...>::SmallIntegerPack(const SmallIntegerPack &other)
    : m_store(other.m_store)
{
}

// ----------------------------------------------------------------------------

template <Uint... Sizes>
SmallIntegerPack<Sizes...> &SmallIntegerPack<Sizes...>::operator=(const SmallIntegerPack &other)
{
  m_store = other.m_store;
  return *this;
}

// ----------------------------------------------------------------------------

template <Uint... Sizes>
SmallIntegerPack<Sizes...>::~SmallIntegerPack()
{
}

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif
