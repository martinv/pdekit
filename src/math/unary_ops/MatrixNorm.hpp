#ifndef PDEKIT_Math_Matrix_Norm_hpp
#define PDEKIT_Math_Matrix_Norm_hpp

#include "math/Matrix.hpp"
#include "math/unary_ops/matrix_unary_ops/MatrixExpressionNorm.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

namespace detail
{

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
struct EuclideanMatrixNormEvaluator
{
  /// Resulting value type
  using value_type = typename MT::value_type;

  /// Type of matrix argument: can be directly the matrix, or in case
  /// the input type MT is in fact an expression (for example matrix-matrix
  /// product), then MatrixArg will be a MT::tensor type, which is typically
  /// DynamicMatrix, or StaticMatrix. A temporary variable of this type will
  /// cache the resulting matrix values (see line 'MatrixArg matrix(mat)' ),
  /// from which the norm will then be computed.
  using MatrixArg = typename common::SelectType<MT::evaluates_fast, const MT &,
                                                const typename MT::tensor_type>::type;

  inline static value_type const evaluate(const MT &mat)
  {
    MatrixArg matrix(mat);

    value_type norm = value_type();

    for (Uint i = 0; i < matrix.rows(); ++i)
    {
      for (Uint j = 0; j < matrix.cols(); ++j)
      {
        norm += matrix(i, j) * matrix(i, j);
      }
    }

    return std::sqrt(norm);
  }
};

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
struct MaxMatrixNormEvaluator
{
  /// Resulting value type
  using value_type = typename MT::value_type;

  /// Type of matrix argument: can be directly the matrix, or in case
  /// the input type MT is in fact an expression (for example matrix-matrix
  /// product), then MatrixArg will be a MT::tensor type, which is typically
  /// DynamicMatrix, or StaticMatrix. A temporary variable of this type will
  /// cache the resulting matrix values (see line 'MatrixArg matrix(mat)' ),
  /// from which the norm will then be computed.
  using MatrixArg = typename common::SelectType<MT::evaluates_fast, const MT &,
                                                const typename MT::tensor_type>::type;

  inline static value_type const evaluate(const MT &mat)
  {
    MatrixArg matrix(mat);

    value_type norm = value_type();

    for (Uint i = 0; i < matrix.rows(); ++i)
    {
      for (Uint j = 0; j < matrix.cols(); ++j)
      {
        norm = std::max(norm, std::abs(matrix(i, j)));
      }
    }

    return norm;
  }
};

// ----------------------------------------------------------------------------

} // namespace detail

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
typename MatrixExpressionNorm<MT, SO, detail::EuclideanMatrixNormEvaluator<MT, SO>>::value_type norm_e2(
    math::DenseMatrix<MT, SO> const &mat)
{
  MatrixExpressionNorm<MT, SO, detail::EuclideanMatrixNormEvaluator<MT, SO>> norm(
      mat.wrapped_type());
  return norm.value();
}

// ----------------------------------------------------------------------------

template <typename MT, bool SO>
typename MatrixExpressionNorm<MT, SO, detail::MaxMatrixNormEvaluator<MT, SO>>::value_type norm_max(
    math::DenseMatrix<MT, SO> const &mat)
{
  MatrixExpressionNorm<MT, SO, detail::MaxMatrixNormEvaluator<MT, SO>> norm(mat.wrapped_type());
  return norm.value();
}

// ----------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif // Matrix_Norm_hpp
