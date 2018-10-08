#ifndef PDEKIT_Math_Scalar_DVec_Mult_hpp
#define PDEKIT_Math_Scalar_DVec_Mult_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"
#include "math/traits/MultTrait.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename ST, typename VT, bool TF>
class DScalarDVecMultExpr : public DenseVector<DScalarDVecMultExpr<ST, VT, TF>, TF>
{
  private:
  /// Result type of the left-hand side scalar expression
  /// This always has to be ScalarConstant< ... >, not an expression
  /// using STT = typename ST::tensor_type;

  /// Result type of the right-hand side dense vector expression
  /// This always has to be Vector< ... >, not an expression
  using VTT = typename VT::tensor_type;

  /// Value type of the left-hand side dense matrix expression
  // using SVT = typename STT::value_type;

  /// Value type of the right-hand side dense vector expression
  using VVT = typename VTT::value_type;

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename MultTrait<ST, VTT>::type;

  /// Data type for composite expression templates
  using composite_type = const DScalarDVecMultExpr &;

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
  /// The scalar - vector product does NOT own data in the sense that the
  /// resulting vector is not explicitly stored.
  enum
  {
    owns_data = 0
  };

  /// Member data type of the left-hand and right-hand side dense vector
  /// expression
  /// Each entry of the rhs vector will be used precisely once in the
  /// scalar-vector multiplication and we would not win anything in
  /// pre-computing and caching the resulting vector in case
  /// VT is not a vector, but an expression resulting into a vector
  using Lhs = const ST; /// The scalar is held by value!

  // using Rhs = const VT&;
  using Rhs = typename common::SelectType<VT::evaluates_fast, VT const &, VTT const>::type;

  /// Constructor
  DScalarDVecMultExpr(const ST &scalar_left, const VT &vec_right)
      : m_scalar(scalar_left), m_vec(vec_right)
  {
  }

  /// Indexing operator
  inline const value_type operator[](const Uint i) const
  {
    return m_scalar * m_vec[i];
  }

  /// Barton-Nackman trick: consider this as overload of the default 'assign'
  /// free function that is defined in DMatEvalExpr.hpp
  template <typename VT2>
  friend inline void free_assign(DenseVector<VT2, TF> &lhs, const DScalarDVecMultExpr &rhs)
  {
    DScalarDVecMultExpr::assign(lhs.wrapped_type(), rhs.m_scalar, rhs.m_vec);
  }

  template <typename VecType1, typename ScalarType, typename VecType2>
  inline static void assign(VecType1 &vec_result, const ScalarType c, const VecType2 &vec_rhs)
  {
    for (Uint i = 0; i < vec_result.size(); ++i)
    {
      vec_result[i] = c * vec_rhs[i];
    }
  }

  /// Length of the resulting vector
  inline Uint size() const
  {
    return m_vec.size();
  }

  private:
  Lhs m_scalar;
  Rhs m_vec;
};

// ----------------------------------------------------------------------------

// template<typename RV, bool TF>
// inline const DScalarDVecProdExpr<ScalarConstant<Real>,RV,TF>
// operator* ( const Real scalar, const DenseVector<RV,TF>& vector )
//{
//  return DScalarDVecProdExpr<ScalarConstant<Real>,RV,TF> (
// ScalarConstant<Real>(scalar),vector.wrapped_type() );
//}

template <typename RV, bool TF>
inline const DScalarDVecMultExpr<Real, RV, TF> operator*(const Real scalar,
                                                         const DenseVector<RV, TF> &vector)
{
  return DScalarDVecMultExpr<Real, RV, TF>(scalar, vector.wrapped_type());
}

// ----------------------------------------------------------------------------
// TENSOR RANK OF VECTOR = 1
// ----------------------------------------------------------------------------

template <typename ST, typename VT, bool TF>
struct TensorRank<DScalarDVecMultExpr<ST, VT, TF>>
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
