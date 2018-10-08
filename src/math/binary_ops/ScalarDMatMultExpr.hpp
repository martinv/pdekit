#ifndef PDEKIT_Math_Scalar_DMat_Mult_hpp
#define PDEKIT_Math_Scalar_DMat_Mult_hpp

#include "common/Meta.hpp"
#include "math/Matrix.hpp"
#include "math/ScalarConstant.hpp"
#include "math/TensorRank.hpp"
#include "math/traits/MultTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename ST, typename MT, bool SO>
class DScalarDMatMultExpr : public DenseMatrix<DScalarDMatMultExpr<ST, MT, SO>, SO>
{
  private:
  // Result type of the left-hand side scalar expression
  // This always has to be ScalarConstant< ... >, not an expression
  // using STT = typename ST::tensor_type;

  /// Result type of the right-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT = typename MT::tensor_type;

  // Value type of the left-hand side dense matrix expression
  // using SVT = typename STT::value_type;

  /// Value type of the right-hand side dense matrix expression
  using MVT = typename MTT::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename MultTrait<ST, MTT>::type;

  /// Data type for composite expression templates
  using composite_type = const DScalarDMatMultExpr &;

  /// Resulting value type
  using value_type = typename tensor_type::value_type;

  enum
  {
    is_expression = 1
  };

  enum
  {
    evaluates_fast = 1
  };
  /// The scalar - matrix product does NOT own data in the sense that the
  /// resulting vector is not explicitly stored.
  enum
  {
    owns_data = 0
  };

  /// Member data type of the left-hand and right-hand side dense matrix
  /// expression
  /// Each entry of the rhs matrix will be used precisely once in the
  /// scalar-matrix multiplication and we would not win anything in
  /// pre-computing and caching the resulting matrix in case MT is not a
  /// vector, but an expression resulting into a matrix
  using Lhs = const ST; /// The scalar is held by value!

  // using Rhs = const MT&;
  using Rhs = typename common::SelectType<MT::evaluates_fast, MT const &, MTT const>::type;

  /// Constructor
  DScalarDMatMultExpr(const ST &scalar_left, const MT &mat_right)
      : m_scalar(scalar_left), m_mat(mat_right)
  {
  }

  /// Indexing operator
  inline const value_type operator()(const Uint i, const Uint j) const
  {
    return m_scalar * m_mat(i, j);
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename MT2>
  friend inline void free_assign(DenseMatrix<MT2, SO> &lhs, const DScalarDMatMultExpr &rhs)
  {
    DScalarDMatMultExpr::assign(lhs.wrapped_type(), rhs.m_scalar, rhs.m_mat);
  }

  template <typename MatType1, typename ScalarType, typename MatType2>
  inline static void assign(MatType1 &mat_result, const ScalarType c, const MatType2 &mat_rhs)
  {
    for (Uint i = 0; i < mat_result.rows(); ++i)
    {
      for (Uint j = 0; j < mat_result.cols(); ++j)
      {
        mat_result(i, j) = c * mat_rhs(i, j);
      }
    }
  }

  /// Number of rows of the resulting matrix
  inline Uint rows() const
  {
    return m_mat.rows();
  }

  /// Number of columns of the resulting matrix
  inline Uint cols() const
  {
    return m_mat.cols();
  }

  private:
  Lhs m_scalar;
  Rhs m_mat;
};

// ----------------------------------------------------------------------------

// template<typename RM, bool SO>
// inline const DScalarDMatProdExpr< ScalarConstant<Real>,RM,SO>
// operator* ( const Real scalar, const DenseMatrix<RM,SO>& matrix )
//{
//  return DScalarDMatProdExpr<ScalarConstant<Real>,RM,SO> (
// ScalarConstant<Real>(scalar),matrix.wrapped_type() );
//}

template <typename RM, bool SO>
inline const DScalarDMatMultExpr<Real, RM, SO> operator*(const Real scalar,
                                                         const DenseMatrix<RM, SO> &matrix)
{
  return DScalarDMatMultExpr<Real, RM, SO>(scalar, matrix.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename ST, typename RM, bool SO>
struct TensorRank<DScalarDMatMultExpr<ST, RM, SO>>
{
  enum
  {
    value = tensor_rank_2
  };
};

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
