#ifndef PDEKIT_Math_Add_Trait_hpp
#define PDEKIT_Math_Add_Trait_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

/// Forward declarations

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct AddTrait;

// ----------------------------------------------------------------------------

// Uint + Uint
template <>
struct AddTrait<Uint, Uint>
{
  using type = Uint;
};

// Uint + Real
template <>
struct AddTrait<Uint, Real>
{
  using type = Real;
};

// Real + Uint
template <>
struct AddTrait<Real, Uint>
{
  using type = Real;
};

// Real + Real
template <>
struct AddTrait<Real, Real>
{
  using type = Real;
};

} // namespace math

} // namespace pdekit

#endif
