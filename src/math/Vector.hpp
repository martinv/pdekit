#ifndef PDEKIT_Math_Vector_hpp
#define PDEKIT_Math_Vector_hpp

namespace pdekit
{

namespace math
{

/*
template<typename Expr>
class ExpressionBase
{
  public:

    using cExpr = const Expr&;

    inline cExpr getCExpr() const
    {
      return static_cast<cExpr> ( *this );
    }

    virtual ~ExpressionBase() { }
};
*/

/// ---------------------------------------------------------------------------
/// Template class to mark a class 'VT' as vector
/// VT ... vector type
/// TF ... transpose flag
/// ---------------------------------------------------------------------------

template <typename VT, bool TF>
class Vector
{
  public:
  /// Typedef for the vector type wrapped by this class
  using VectorType = VT;

  /// Return this converted to wrapped type
  inline VectorType &wrapped_type()
  {
    return static_cast<VectorType &>(*this);
  }

  /// Return this instance converted to wrapped type, const version
  inline const VectorType &wrapped_type() const
  {
    return static_cast<const VectorType &>(*this);
  }
};

/// ---------------------------------------------------------------------------
/// Template class to mark a class 'VT' as dense vector
/// VT ... vector type
/// TF ... transpose flag
/// Currently, a dense vector can be the class StaticVector, DynamicVector,
/// VectorBlock or any other class that represents a mathematical expression
/// whose result is a dense vector
/// (for example, a sum DynamicVector + DynamicVector, or StaticVector +
/// VectorBlock
/// or the product constant * VectorBlock )
/// ---------------------------------------------------------------------------

template <typename VT, bool TF>
struct DenseVector : public Vector<VT, TF>
{
};

// ---------------------------------------------------------------------------

} // namespace math

} // namespace pdekit

#endif
