#ifndef PDEKIT_Math_DMat_DMat_Add_hpp
#define PDEKIT_Math_DMat_DMat_Add_hpp

#include "common/Meta.hpp"
#include "math/Matrix.hpp"
#include "math/TensorRank.hpp"
#include "math/traits/AddTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT1, typename MT2, bool SO>
class DMatDMatAddExpr : public DenseMatrix<DMatDMatAddExpr<MT1, MT2, SO>, SO>
{
  private:
  /// Result type of the left-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT1 = typename MT1::tensor_type;

  /// Result type of the right-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT2 = typename MT2::tensor_type;

  /// Value type of the left-hand side dense matrix expression
  using MVT1 = typename MTT1::value_type;

  /// Value type of the right-hand side dense matrix expression
  using MVT2 = typename MTT2::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename AddTrait<MTT1, MTT2>::type;

  /// Data type for composite expression templates
  using composite_type = const DMatDMatAddExpr &;

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
  enum
  {
    owns_data = 0
  };

  /// Member data type of the left-hand side dense matrix expression
  /// In case the Lhs is a result of some more complex expression which is
  /// slow to evaluate, Lhs will be of type 'DynamicMatrix/StaticMatrix' and
  /// cache the result of this expression. If the lhs operand is itself a
  /// vector or is fast to evaluate (e.g. addition of two vectors), then Lhs
  /// type will be of type 'reference to operand' Similarly for Rhs
  using Lhs = typename common::SelectType<MT1::evaluates_fast, MT1 const &, MTT1 const>::type;

  using Rhs = typename common::SelectType<MT2::evaluates_fast, MT2 const &, MTT2 const>::type;

  /// Constructor
  DMatDMatAddExpr(const MT1 &lhs, const MT2 &rhs) : m_lhs(lhs), m_rhs(rhs)
  {
  }

  /// Operator returning one value
  inline const value_type operator()(const Uint i, const Uint j) const
  {
    return m_lhs(i, j) + m_rhs(i, j);
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename MT3>
  friend inline void free_assign(DenseMatrix<MT3, SO> &lhs, const DMatDMatAddExpr &rhs)
  {
    DMatDMatAddExpr::assign(lhs.wrapped_type(), rhs.m_lhs, rhs.m_rhs);
  }

  template <typename MatType1, typename MatType2, typename MatType3>
  inline static void assign(MatType1 &mat_result, const MatType2 &mat_left,
                            const MatType3 &mat_right)
  {
    for (Uint i = 0; i < mat_result.rows(); ++i)
    {
      for (Uint j = 0; j < mat_result.cols(); ++j)
      {
        mat_result(i, j) = mat_left(i, j) + mat_right(i, j);
      }
    }
  }

  /// Return the number of rows of the resulting matrix
  inline Uint rows() const
  {
    return m_lhs.rows();
  }

  /// Return the number of columns of the resulting matrix
  inline Uint cols() const
  {
    return m_lhs.cols();
  }

  private:
  Lhs m_lhs;
  Rhs m_rhs;
};

// ----------------------------------------------------------------------------

/// Operator for addition of two dense matrices
template <typename LM, typename RM, bool SO>
inline const DMatDMatAddExpr<LM, RM, SO> operator+(const DenseMatrix<LM, SO> &lhs,
                                                   const DenseMatrix<RM, SO> &rhs)
{
  return DMatDMatAddExpr<LM, RM, SO>(lhs.wrapped_type(), rhs.wrapped_type());
}

// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// TENSOR RANK OF MATRIX = 2
// ----------------------------------------------------------------------------

template <typename LM, typename RM, bool SO>
struct TensorRank<DMatDMatAddExpr<LM, RM, SO>>
{
  enum
  {
    value = tensor_rank_2
  };
};

} // Namespace math

} // Namespace pdekit

#endif
