#ifndef PDEKIT_Math_Math_Op_Traits_hpp
#define PDEKIT_Math_Math_Op_Traits_hpp

#include "common/PDEKit.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename T1, typename T2>
struct MathOpTrait;

// ----------------------------------------------------------------------------

/// Binary operation OP(Uint,Uint) traits:
template <>
struct MathOpTrait<Uint, Uint>
{
  // Common high data type
  using high_type = Uint;

  // Common low data type
  using low_type = Uint;
};

// ----------------------------------------------------------------------------

/// Binary operation OP(Uint,Real) traits:
template <>
struct MathOpTrait<Uint, Real>
{
  // Common high data type
  using high_type = Real;

  // Common low data type
  using low_type = Uint;
};

// ----------------------------------------------------------------------------

/// Binary operation OP(Real,Uint) traits:
template <>
struct MathOpTrait<Real, Uint>
{
  // Common high data type
  using high_type = Real;

  // Common low data type
  using low_type = Uint;
};

// ----------------------------------------------------------------------------

/// Binary operation OP(Real,Real) traits:
template <>
struct MathOpTrait<Real, Real>
{
  // Common high data type
  using high_type = Real;

  // Common low data type
  using low_type = Real;
};

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
