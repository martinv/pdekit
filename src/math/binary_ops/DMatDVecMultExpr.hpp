#ifndef PDEKIT_Math_DMat_DVec_Mult_hpp
#define PDEKIT_Math_DMat_DVec_Mult_hpp

#include "common/Meta.hpp"
#include "math/Matrix.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"
#include "math/traits/MultTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT, typename VT, bool SO, bool TF>
class DMatDVecMultExpr : public DenseVector<DMatDVecMultExpr<MT, VT, SO, TF>, TF>
{
  private:
  /// Result type of the left-hand side dense matrix expression
  /// This always has to be Matrix< ... >, not an expression
  using MTT = typename MT::tensor_type;

  /// Result type of the right-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT = typename VT::tensor_type;

  /// Value type of the left-hand side dense matrix expression
  using MVT = typename MTT::value_type;

  /// Value type of the right-hand side dense vector expression
  using VVT = typename VTT::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename MultTrait<MTT, VTT>::type;

  /// Data type for composite expression templates
  using composite_type = const DMatDVecMultExpr &;

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

  /// The matrix - vector product does NOT own data in the sense that the
  /// resulting vector is not explicitly stored.
  enum
  {
    owns_data = 0
  };

  /// Member data type of the left-hand side dense matrix expression
  /// We don't need to decide whether this should be (matrix) value or
  /// reference to matrix Each entry of the matrix will be used precisely once
  /// in the matrix-vector multiplication and we would not win anything in
  /// pre-computing and caching the matrix in case MT is not a matrix but a
  /// matrix expression (e.g. sum or product of two matrices)
  using Lhs = typename MT::composite_type;

  // using Rhs = typename common::SelectType<VT::is_expression, const VTT,
  // const VT&>::type;
  using Rhs = typename common::SelectType<VT::evaluates_fast, const VT &, const VTT>::type;

  /// Constructor
  DMatDVecMultExpr(const MT &mat, const VT &vec) : m_mat(mat), m_vec(vec), N(m_mat.cols())
  {
  }

  /// Element access operator
  inline const value_type operator[](const Uint i) const
  {
    // return
    // vector_vector_product_impl<value_type,Lhs,CTV,MT::is_expression,VT::is_expression>::apply
    // ( m_mat,m_vec,N,i*N,0 );

    value_type result = m_mat(i, 0) * m_vec[0];
    for (Uint j = 1; j < m_mat.cols(); ++j)
    {
      result += m_mat(i, j) * m_vec[j];
    }
    return result;
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename VT2>
  friend inline void free_assign(DenseVector<VT2, SO> &lhs, const DMatDVecMultExpr &rhs)
  {
    DMatDVecMultExpr::assign(lhs.wrapped_type(), rhs.m_mat, rhs.m_vec);
  }

  template <typename VT1, typename VT2, typename MT1>
  inline static void assign(VT1 &vec_result, const MT1 &mat_left, const VT2 &vec_right)
  {
    for (Uint i = 0; i < mat_left.rows(); ++i)
    {
      vec_result[i] = mat_left(i, 0) * vec_right[0];
      for (Uint j = 1; j < mat_left.cols(); ++j)
      {
        vec_result[i] += mat_left(i, j) * vec_right[j];
      }
    }
  }

  /// Length of the resulting vector
  inline Uint size() const
  {
    return m_mat.rows();
  }

  private:
  /// Composite type of the vector expression: if the rhs is a vector, then
  /// the CTV will be reference to a sparse/dense vector or vectorview
  /// (depending on what type is VTT), otherwise CTV will be reference to VT -
  /// the vector expression
  using CTV = typename common::SelectType<VT::is_expression, const VTT &, const VT &>::type;

  /// Operands
  Lhs m_mat;
  Rhs m_vec;

  /// Length of the resulting vector
  const Uint N;
};

// ----------------------------------------------------------------------------

/// Operator for dense matrix - dense vector multiplication:
template <typename LM, typename RV, bool SO, bool TF>
inline const DMatDVecMultExpr<LM, RV, SO, TF> operator*(const DenseMatrix<LM, SO> &mat,
                                                        const DenseVector<RV, TF> &vec)
{
  return DMatDVecMultExpr<LM, RV, SO, TF>(mat.wrapped_type(), vec.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename LM, typename RV, bool SO, bool TF>
struct TensorRank<DMatDVecMultExpr<LM, RV, SO, TF>>
{
  enum
  {
    value = tensor_rank_1
  };
};

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
