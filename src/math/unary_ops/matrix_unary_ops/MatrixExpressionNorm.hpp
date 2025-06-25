#ifndef PDEKIT_Math_Matrix_Norm_Expression_hpp
#define PDEKIT_Math_Matrix_Norm_Expression_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT, bool SO, typename NormEvaluator>
class MatrixExpressionNorm
{

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename NormEvaluator::value_type;

  /// Data type for composite expression templates
  using composite_type = const MatrixExpressionNorm &;

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
  MatrixExpressionNorm(const MT &mat)
  {
    m_norm = NormEvaluator::evaluate(mat);
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

template <typename MT, bool SO, typename NormEvaluator>
struct TensorRank<MatrixExpressionNorm<MT, SO, NormEvaluator>>
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
