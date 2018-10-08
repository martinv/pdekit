#ifndef PDEKIT_Math_DVecDVecAdd_hpp
#define PDEKIT_Math_DVecDVecAdd_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"
#include "math/traits/AddTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename VT1, typename VT2, bool TF>
class DVecDVecAddExpr : public DenseVector<DVecDVecAddExpr<VT1, VT2, TF>, TF>
{
  private:
  /// Result type of the left-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT1 = typename VT1::tensor_type;

  /// Result type of the right-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT2 = typename VT2::tensor_type;

  /// Value type of the left-hand side dense vector expression
  using VVT1 = typename VTT1::value_type;

  /// Value type of the right-hand side dense vector expression
  using VVT2 = typename VTT2::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename AddTrait<VTT1, VTT2>::type;

  /// Data type for composite expression templates
  using composite_type = const DVecDVecAddExpr &;

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
  /// The matrix addition does NOT own data in the sense that the resulting
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

  // Operand type selection logic:
  // If the Lhs evaluates fast,
  //                      then the Lhs type will be 'reference to const' in
  //                      case it owns data or it is an expression If it does
  //                      not own data and it is not an expression, it must be
  //                      a proxy object such as BlockVector or
  //                      ConstBlockVector which is stored by value (this is
  //                      useful in interpolation::VectorMeshFunction and its
  //                      algebraic operators where 'entry_types' are
  //                      VectorBlocks returned by methods as temporary
  //                      variables
  // else if the Lhs evaluates slow, then the Lhs is a static/dynamic vector
  // and the intermediate result will be cached.

  /*
  using Lhs = typename common::SelectType<
      VT1::evaluates_fast, typename common::SelectType < VT1::owns_data ||
  VT1::is_expression, VT1 const&, VT1 const>::type, VTT1 const > ::type;

  using Rhs = typename common::SelectType<
      VT2::evaluates_fast, typename common::SelectType < VT2::owns_data ||
  VT2::is_expression, VT2 const&, VT2 const>::type, VTT2 const > ::type;
  */

  /// Constructor
  DVecDVecAddExpr(const VT1 &vec_left, const VT2 &vec_right) : m_vec1(vec_left), m_vec2(vec_right)
  {
  }

  /// Operator returning one value
  inline const value_type operator[](const Uint i) const
  {
    return m_vec1[i] + m_vec2[i];
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename VT3>
  friend inline void free_assign(DenseVector<VT3, TF> &lhs, const DVecDVecAddExpr &rhs)
  {
    DVecDVecAddExpr::assign(lhs.wrapped_type(), rhs.m_vec1, rhs.m_vec2);
  }

  template <typename VecType1, typename VecType2, typename VecType3>
  inline static void assign(VecType1 &vec_result, const VecType2 &vec_left,
                            const VecType3 &vec_right)
  {
    for (Uint i = 0; i < vec_result.size(); ++i)
    {
      vec_result[i] = vec_left[i] + vec_right[i];
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

/// Operator for addition of two dense vectors:
template <typename LV, typename RV, bool TF>
inline const DVecDVecAddExpr<LV, RV, TF> operator+(const DenseVector<LV, TF> &lhs,
                                                   const DenseVector<RV, TF> &rhs)
{
  return DVecDVecAddExpr<LV, RV, TF>(lhs.wrapped_type(), rhs.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename LV, typename RV, bool TF>
struct TensorRank<DVecDVecAddExpr<LV, RV, TF>>
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
