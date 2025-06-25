#ifndef PDEKIT_Math_Vector_Norm_Expression_hpp
#define PDEKIT_Math_Vector_Norm_Expression_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename VT, bool TF, typename NormEvaluator>
class VectorNormExpression
{

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename NormEvaluator::value_type;

  /// Data type for composite expression templates
  using composite_type = const VectorNormExpression &;

  /// Resulting value type
  using value_type = typename NormEvaluator::value_type;

  enum
  {
    is_expression = 1
  };

  enum
  {
    evaluates_fast = 1
  };

  /// Constructor
  VectorNormExpression(const VT &vec)
  {
    m_norm = NormEvaluator::evaluate(vec);
  }

  /// Element access operator
  inline const value_type value() const
  {
    return m_norm;
  }

  /// Length of the resulting expression
  inline Uint size() const
  {
    return 1u;
  }

  private:
  /// value of the norm
  value_type m_norm;
};

// ----------------------------------------------------------------------------
// TENSOR RANK OF NORM(SCALAR) = 0
// ----------------------------------------------------------------------------

template <typename VT, bool TF, typename NormEvaluator>
struct TensorRank<VectorNormExpression<VT, TF, NormEvaluator>>
{
  enum
  {
    value = tensor_rank_0
  };
};

// ----------------------------------------------------------------------------

} // Namespace math

} // Namespace pdekit

#endif
