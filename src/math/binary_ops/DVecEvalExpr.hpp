#ifndef PDEKIT_Math_DVecEvalExpr_hpp
#define PDEKIT_Math_DVecEvalExpr_hpp

#include "math/Vector.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename VT1, typename VT2, bool SO>
inline void free_assign(DenseVector<VT1, SO> &v_lhs, const DenseVector<VT2, SO> &v_rhs)
{
  // TODO: check that the number of rows of lhs and rhs is equal
  v_lhs.wrapped_type().assign(v_rhs);
}

// ----------------------------------------------------------------------------
} // namespace math
} // namespace pdekit

#endif
