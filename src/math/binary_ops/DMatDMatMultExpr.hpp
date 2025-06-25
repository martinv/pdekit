#ifndef PDEKIT_Math_DMat_DMat_Mult_hpp
#define PDEKIT_Math_DMat_DMat_Mult_hpp

#include "common/Meta.hpp"
#include "math/Matrix.hpp"
#include "math/TensorRank.hpp"
#include "math/traits/MultTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT1, typename MT2, bool SO>
class DMatDMatMultExpr : public DenseMatrix<DMatDMatMultExpr<MT1, MT2, SO>, SO>
{
  private:
  /// Result type of the left-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT1 = typename MT1::tensor_type;

  /// Result type of the right-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT2 = typename MT2::tensor_type;

  /// Composite type for the left-hand side dense matrix expression
  using CT1 = typename MT1::composite_type;

  /// Composite type of the right-hand side dense matrix expression
  using CT2 = typename MT2::composite_type;

  /// Value type of the left-hand side dense matrix expression
  using MVT1 = typename MTT1::value_type;

  /// Value type of the right-hand side dense matrix expression
  using MVT2 = typename MTT2::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename MultTrait<MTT1, MTT2>::type;

  /// Data type for composite expression templates
  using composite_type = const DMatDMatMultExpr &;

  /// Resulting value type
  using value_type = typename tensor_type::value_type;

  enum
  {
    is_expression = 1
  };

  enum
  {
    evaluates_fast = 0
  };

  /// The matrix product does NOT own data in the sense that the resulting
  /// matrix is not explicitly stored.
  enum
  {
    owns_data = 0
  };

  /// If the left operand is an expression, then the member variable type Lhs
  /// will be dense matrix, which will cache the result of this expression.
  /// Otherwise (left operand itself is a dense matrix, not matrix
  /// expresssion), the Lhs type will be the 'composite type' - a reference to
  /// the left operand (matrix)
  using Lhs = typename common::SelectType<MT1::is_expression, const MTT1, CT1>::type;

  using Rhs = typename common::SelectType<MT2::is_expression, const MTT2, CT2>::type;

  /// Constructor
  DMatDMatMultExpr(const MT1 &lhs, const MT2 &rhs) : m_lhs(lhs), m_rhs(rhs) //,
  // M(mat_left.rows()),
  // N(mat_left.cols()),
  // P(mat_right.cols())
  {
  }

  // ----------------------------------------------------------------------------

  inline const value_type operator()(const Uint i, const Uint j) const
  {
    value_type result = m_lhs(i, 0) * m_rhs(0, j);
    for (Uint k = 1; k < m_lhs.cols(); ++k)
    {
      result += m_lhs(i, k) * m_rhs(k, j);
    }
    return result;
  }

  // ----------------------------------------------------------------------------

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename MT3>
  friend inline void free_assign(DenseMatrix<MT3, SO> &lhs, const DMatDMatMultExpr &rhs)
  {
    DMatDMatMultExpr::assign(lhs.wrapped_type(), rhs.m_lhs, rhs.m_rhs);
  }

  template <typename MatType1, typename MatType2, typename MatType3>
  inline static void assign(MatType1 &mat_result, const MatType2 &mat_left,
                            const MatType3 &mat_right)
  {
    for (Uint i = 0; i < mat_left.rows(); ++i)
    {
      for (Uint j = 0; j < mat_right.cols(); ++j)
      {
        mat_result(i, j) = mat_left(i, 0) * mat_right(0, j);
      }

      for (Uint k = 1; k < mat_left.cols(); ++k)
      {
        for (Uint j = 0; j < mat_right.cols(); ++j)
        {
          // C(i,j) += A(i,k) * B(k,j)
          // ptr = (target_ptr+i*m_rhs.cols()+j);
          mat_result(i, j) += mat_left(i, k) * mat_right(k, j);
        }
      }
    } // Loop over rows of lhs
  }

  // ----------------------------------------------------------------------------

  inline Uint size() const
  {
    return m_lhs.rows() * m_rhs.cols(); // M*P;
  }

  // ----------------------------------------------------------------------------

  inline Uint rows() const
  {
    return m_lhs.rows(); // M;
  }

  // ----------------------------------------------------------------------------

  inline Uint cols() const
  {
    return m_rhs.cols(); // P;
  }

  private:
  // Operands: the left-hand side matrix and the
  // right-hand side matrix in the product M1 x M2

  Lhs m_lhs;
  Rhs m_rhs;

  // Sizes of matrices: (M x N) x ( N x P )
  // The size of the resulting matrix product will be M x P
  //    const Uint M;
  //    const Uint N;
  //    const Uint P;
};

// ----------------------------------------------------------------------------

/// Operator for multiplication of two dense matrices:
template <typename LM, typename RM, bool SO>
inline const DMatDMatMultExpr<LM, RM, SO> operator*(const DenseMatrix<LM, SO> &lhs_mat,
                                                    const DenseMatrix<RM, SO> &rhs_mat)
{
  return DMatDMatMultExpr<LM, RM, SO>(lhs_mat.wrapped_type(), rhs_mat.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename LM, typename RM, bool SO>
struct TensorRank<DMatDMatMultExpr<LM, RM, SO>>
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
