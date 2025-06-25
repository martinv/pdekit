#ifndef PDEKIT_Math_Vector_Norm_hpp
#define PDEKIT_Math_Vector_Norm_hpp

#include "math/Vector.hpp"
#include "math/unary_ops/vector_unary_ops/VectorNormExpression.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

namespace detail
{

// ----------------------------------------------------------------------------

template <typename VT, bool TF>
struct EuclideanVectorNormEvaluator
{
  /// Resulting value type
  using value_type = typename VT::value_type;

  /// Type of vector argument: can be directly the vector, or in case
  /// the input type VT is in fact an expression (for example matrix-vector
  /// product), then VectorArg will be a VT::tensor type, which is typically
  /// DynamicVector, or StaticVector. A temporary variable of this type will
  /// cache the resulting vector values (see line 'VectorArg vector(vec)' ),
  /// from which the norm will then be computed.
  using VectorArg = typename common::SelectType<VT::evaluates_fast, const VT &,
                                                const typename VT::tensor_type>::type;

  inline static value_type const evaluate(const VT &vec)
  {
    VectorArg vector(vec);

    value_type norm = value_type();

    for (Uint i = 0; i < vector.size(); ++i)
    {
      norm += vector[i] * vector[i];
    }

    return std::sqrt(norm);
  }
};

// ----------------------------------------------------------------------------

template <typename VT, bool TF>
struct MaxVectorNormEvaluator
{
  /// Resulting value type
  using value_type = typename VT::value_type;

  /// Type of vector argument: can be directly the vector, or in case
  /// the input type VT is in fact an expression (for example matrix-vector
  /// product), then VectorArg will be a VT::tensor type, which is typically
  /// DynamicVector, or StaticVector. A temporary variable of this type will
  /// cache the resulting vector values (see line 'VectorArg vector(vec)' ),
  /// from which the norm will then be computed.
  using VectorArg = typename common::SelectType<VT::evaluates_fast, const VT &,
                                                const typename VT::tensor_type>::type;

  inline static value_type const evaluate(const VT &vec)
  {
    VectorArg vector(vec);

    value_type norm = value_type();

    for (Uint i = 0; i < vector.size(); ++i)
    {
      norm = std::max(norm, std::abs(vector[i]));
    }

    return norm;
  }
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

template <typename VT, bool TF>
typename VectorNormExpression<VT, TF, detail::EuclideanVectorNormEvaluator<VT, TF>>::value_type norm_e2(
    math::DenseVector<VT, TF> const &vec)
{
  VectorNormExpression<VT, TF, detail::EuclideanVectorNormEvaluator<VT, TF>> norm(
      vec.wrapped_type());
  return norm.value();
}

// ----------------------------------------------------------------------------

template <typename VT, bool TF>
typename VectorNormExpression<VT, TF, detail::MaxVectorNormEvaluator<VT, TF>>::value_type norm_max(
    math::DenseVector<VT, TF> const &vec)
{
  VectorNormExpression<VT, TF, detail::MaxVectorNormEvaluator<VT, TF>> norm(vec.wrapped_type());
  return norm.value();
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif // Vector_Norm_hpp
