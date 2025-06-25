#ifndef PDEKIT_Math_DMatEvalExpr_hpp
#define PDEKIT_Math_DMatEvalExpr_hpp

#include "math/Matrix.hpp"

namespace pdekit
{

namespace math
{

// ----------------------------------------------------------------------------

template <typename MT1, typename MT2, bool SO>
inline void free_assign(DenseMatrix<MT1, SO> &m_lhs, const DenseMatrix<MT2, SO> &m_rhs)
{
  // TODO: check that the number of columns and rows of lhs and rhs is equal
  m_lhs.wrapped_type().assign(m_rhs);
}

// ----------------------------------------------------------------------------
} // namespace math
} // namespace pdekit

#endif
