#ifndef PDEKIT_Math_Mult_Trait_hpp
#define PDEKIT_Math_Mult_Trait_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct MultTrait;

// ----------------------------------------------------------------------------

// Uint * Uint
template <>
struct MultTrait<Uint, Uint>
{
  using type = Uint;
};

// Uint * Real
template <>
struct MultTrait<Uint, Real>
{
  using type = Real;
};

// Real * Uint
template <>
struct MultTrait<Real, Uint>
{
  using type = Real;
};

// Real * Real
template <>
struct MultTrait<Real, Real>
{
  using type = Real;
};

} // namespace math

} // namespace pdekit

#endif
