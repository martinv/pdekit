#ifndef PDEKIT_Math_Div_Trait_hpp
#define PDEKIT_Math_Div_Trait_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct DivTrait;

// ----------------------------------------------------------------------------

// Uint / Uint
template <>
struct DivTrait<Uint, Uint>
{
  using type = Uint;
};

// Uint / Real
template <>
struct DivTrait<Uint, Real>
{
  using type = Real;
};

// Real / Uint
template <>
struct DivTrait<Real, Uint>
{
  using type = Real;
};

// Real / Real
template <>
struct DivTrait<Real, Real>
{
  using type = Real;
};

} // namespace math

} // namespace pdekit

#endif
