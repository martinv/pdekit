#ifndef PDEKIT_Common_Meta_hpp
#define PDEKIT_Common_Meta_hpp

#include <string>

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace common
{

// ----------------------------------------------------------------------------

struct NullType
{
};

struct TrueType
{
};

struct FalseType
{
};

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct TypesAreIdentical
{
  enum
  {
    value = 0
  };
};

template <typename T>
struct TypesAreIdentical<T, T>
{
  enum
  {
    value = 1
  };
};

// ----------------------------------------------------------------------------

template <bool C, typename T1, typename T2>
struct SelectType
{
  typedef T1 type;
};

template <typename T1, typename T2>
struct SelectType<false, T1, T2>
{
  typedef T2 type;
};

// ----------------------------------------------------------------------------

template <bool C, int I1, int I2>
struct SelectValue
{
  enum
  {
    value = I1
  };
};

template <int I1, int I2>
struct SelectValue<false, I1, I2>
{
  enum
  {
    value = I2
  };
};

// ----------------------------------------------------------------------------

template <Uint N1, Uint N2>
struct StaticMax
{
  static const Uint value = (N1 > N2) ? N1 : N2;
};

// ----------------------------------------------------------------------------

template <Uint N1, Uint N2>
struct StaticMin
{
  static const Uint value = (N1 > N2) ? N2 : N1;
};

// ----------------------------------------------------------------------------

template <Uint N1, Uint N2>
struct IntegersAreEqual
{
  enum
  {
    value = 0
  };
};

// ----------------------------------------------------------------------------

template <Uint N>
struct IntegersAreEqual<N, N>
{
  enum
  {
    value = 1
  };
};

// ----------------------------------------------------------------------------

template <Uint N>
struct IsEvenInteger
{
  enum
  {
    value = (N % 2 == 0)
  };
};

// ----------------------------------------------------------------------------

// Minimum number of bits that is needed to store integer N
// Note that for example for N = 16, we need 5 bits, because 4 bits
// are only sufficient to store values 0,...,15

template <Uint N>
struct MinNbBitsToStoreNumber
{
  enum
  {
    value = 1 + MinNbBitsToStoreNumber<N / 2>::value
  };
};

// ----------------------------------------------------------------------------

template <>
struct MinNbBitsToStoreNumber<0U>
{
  enum
  {
    value = 1
  };
};

// ----------------------------------------------------------------------------

template <>
struct MinNbBitsToStoreNumber<1U>
{
  enum
  {
    value = 1
  };
};

// ----------------------------------------------------------------------------

// This class holds the value of 2^N
template <Uint N>
struct PowerOfTwo
{
  enum
  {
    value = 2 * PowerOfTwo<N - 1>::value
  };
};

// ----------------------------------------------------------------------------

// This class holds the value of 2^0
template <>
struct PowerOfTwo<0U>
{
  enum
  {
    value = 1
  };
};

// ----------------------------------------------------------------------------

// This class holds a number which is given by
// 0 0 0 0 0 0 1 1 1 1 1 in binary representation.
// The parameter N determines the number of least-significant bits
// that are switched to 1

template <Uint N>
struct BitsToDecimal
{
  enum
  {
    value = PowerOfTwo<N>::value - 1
  };
};

// ----------------------------------------------------------------------------

template <Uint N>
struct PowerOfTwoLargerOrEqualTo
{
  //   enum { value = 2 * PowerOfTwoLargerOrEqualTo<N/2>::value }; // This is
  // not correct the number is twice smaller than it should be
  enum
  {
    value = PowerOfTwo<MinNbBitsToStoreNumber<N>::value>::value
  };
};

// ----------------------------------------------------------------------------

// Sum of a pack of integers. Uses variable template parameters
// Example usage: const unsigned my_sum = Sum<18,8,3,4>::value;

// Forward declaration!
template <Uint...>
struct Sum;

template <Uint size>
struct Sum<size>
{
  enum
  {
    value = size
  };
};

template <Uint size, Uint... sizes>
struct Sum<size, sizes...>
{
  enum
  {
    value = size + Sum<sizes...>::value
  };
};

// ----------------------------------------------------------------------------

} // Namespace common

} // Namespace pdekit

#endif
