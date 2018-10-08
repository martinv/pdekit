#ifndef PDEKIT_Math_Matrix_Expression_Cond_Number_hpp
#define PDEKIT_Math_Matrix_Expression_Cond_Number_hpp

#include "common/Meta.hpp"
#include "math/TensorRank.hpp"
#include "math/Vector.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT, bool SO, typename CondNumberEvaluator>
class MatrixExpressionCondNumber
{

  public:
  /// Result type for expression template evaluation
  using tensor_type = typename CondNumberEvaluator::value_type;

  /// Data type for composite expression templates
  using composite_type = const MatrixExpressionCondNumber &;

  /// Resulting value type
  using value_type = typename CondNumberEvaluator::value_type;

  enum
  {
    is_expression = 1
  };

  enum
  {
    evaluates_fast = 1
  };

  /// Constructor
  MatrixExpressionCondNumber(const MT &mat)
  {
    m_cond_nr = CondNumberEvaluator::evaluate(mat);
  }

  /// Element access operator
  inline const value_type value() const
  {
    return m_cond_nr;
  }

  /// Length of the resulting expression
  inline Uint size() const
  {
    return 1u;
  }

  private:
  /// value of the norm
  value_type m_cond_nr;
};

// ----------------------------------------------------------------------------
// TENSOR RANK OF CONDITION NUMBER (SCALAR) = 0
// ----------------------------------------------------------------------------

template <typename MT, bool SO, typename CondNumberEvaluator>
struct TensorRank<MatrixExpressionCondNumber<MT, SO, CondNumberEvaluator>>
{
  enum
  {
    value = tensor_rank_0
  };
};

// ============================================================================

} // Namespace math

} // Namespace pdekit

#endif
