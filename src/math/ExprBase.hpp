#ifndef PDEKIT_Math_Expr_Base_hpp
#define PDEKIT_Math_Expr_Base_hpp

namespace pdekit
{

namespace math
{

// ============================================================================

// The purpose of the class Scalar is to mark a class as a scalar entity.
// Currently, a scalar vector can be the class ScalarConstant or ScalarVariable

template <typename T>
class Scalar
{
  public:
  inline const T &scalar_expr() const
  {
    return static_cast<const T &>(*this);
  }
};

} // Namespace math

} // Namespace pdekit

#endif
