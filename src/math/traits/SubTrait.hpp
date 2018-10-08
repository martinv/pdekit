#ifndef PDEKIT_Math_Sub_Trait_hpp
#define PDEKIT_Math_Sub_Trait_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct SubTrait;

// ----------------------------------------------------------------------------

// Uint - Uint
template <>
struct SubTrait<Uint, Uint>
{
  using type = Uint;
};

// Uint - Real
template <>
struct SubTrait<Uint, Real>
{
  using type = Real;
};

// Real - Uint
template <>
struct SubTrait<Real, Uint>
{
  using type = Real;
};

// Real - Real
template <>
struct SubTrait<Real, Real>
{
  using type = Real;
};

} // namespace math

} // namespace pdekit

#endif
