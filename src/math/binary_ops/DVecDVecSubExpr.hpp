#ifndef PDEKIT_Math_DVecDVecSub_hpp
#define PDEKIT_Math_DVecDVecSub_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"
#include "math/traits/SubTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename VT1, typename VT2, bool TF>
class DVecDVecSubExpr : public DenseVector<DVecDVecSubExpr<VT1, VT2, TF>, TF>
{
  private:
  /// Result type of the left-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT1 = typename VT1::tensor_type;

  /// Result type of the right-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT2 = typename VT2::tensor_type;

  /// Value type of the left-hand side dense matrix expression
  using VVT1 = typename VTT1::value_type;

  /// Value type of the right-hand side dense vector expression
  using VVT2 = typename VTT2::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename SubTrait<VTT1, VTT2>::type;

  /// Data type for composite expression templates
  using composite_type = const DVecDVecSubExpr &;

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
  /// The vector subtraction does NOT own data in the sense that the resulting
  /// vector is not explicitly stored.
  enum
  {
    owns_data = 0
  };

  /// Member data type of the left-hand side dense vector expression
  /// In case the Lhs is a result of some more complex expression which is
  /// slow to evaluate, Lhs will be of type 'vector' and cache the result of
  /// this expression. If the lhs operand is itself a vector or is fast to
  /// evaluate (e.g. addition of two vectors), then Lhs type will be of type
  /// 'reference to operand' Similarly for Rhs
  using Lhs = typename common::SelectType<VT1::evaluates_fast, VT1 const &, VTT1 const>::type;

  using Rhs = typename common::SelectType<VT2::evaluates_fast, VT2 const &, VTT2 const>::type;

  /// Constructor
  DVecDVecSubExpr(const VT1 &vec_left, const VT2 &vec_right) : m_vec1(vec_left), m_vec2(vec_right)
  {
  }

  /// Operator returning one value
  inline const value_type operator[](const Uint i) const
  {
    return m_vec1[i] - m_vec2[i];
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename VT3>
  friend inline void free_assign(DenseVector<VT3, TF> &lhs, const DVecDVecSubExpr &rhs)
  {
    DVecDVecSubExpr::assign(lhs.wrapped_type(), rhs.m_vec1, rhs.m_vec2);
  }

  template <typename VecType1, typename VecType2, typename VecType3>
  inline static void assign(VecType1 &vec_result, const VecType2 &vec_left,
                            const VecType3 &vec_right)
  {
    for (Uint i = 0; i < vec_result.size(); ++i)
    {
      vec_result[i] = vec_left[i] - vec_right[i];
    }
  }

  /// Return the length of the resulting vector
  inline Uint size() const
  {
    return m_vec1.size();
  }

  private:
  Lhs m_vec1;
  Rhs m_vec2;
};

// ----------------------------------------------------------------------------

/// Operator for subtraction of two dense vectors:
template <typename LV, typename RV, bool TF>
inline const DVecDVecSubExpr<LV, RV, TF> operator-(const DenseVector<LV, TF> &lhs,
                                                   const DenseVector<RV, TF> &rhs)
{
  return DVecDVecSubExpr<LV, RV, TF>(lhs.wrapped_type(), rhs.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename LV, typename RV, bool TF>
struct TensorRank<DVecDVecSubExpr<LV, RV, TF>>
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