#ifndef PDEKIT_Math_Product_Implementation_hpp
#define PDEKIT_Math_Product_Implementation_hpp

namespace pdekit
{

namespace math
{

#if 0

/// vector_vector_product_impl implements a scalar product of 2 vectors
/// - It can be a product row of a matrix * vector, for example
/// generic version:

template<typename ValueType,typename LExprRef, typename RExprRef, Uint LhsIsExpression, Uint RhsIsExpression>
class vector_vector_product_impl
{
  public:
    static inline ValueType apply ( const LExprRef lexpr, const RExprRef rexpr,
                                    const Uint size, const Uint offset_l, const Uint offset_r )
    {
      ValueType result = ValueType();
      for ( Uint i = offset_l, j = offset_r; i < size+offset_l; ++i, ++j )
        result += lexpr[i]*rexpr[j];
      return result;
    }
};

// Partial specialization in case the left operand has direct access to the data - should run a bit faster
// In other words, the left operand is not an expression, while the right one is
// and cannot access the stored data directly, but only through the [] operator

template<typename ValueType,typename LExprRef, typename RExprRef>
class vector_vector_product_impl<ValueType,LExprRef,RExprRef,0,1>
{
  public:
    static inline ValueType apply ( const LExprRef lexpr, const RExprRef rexpr,
                                    const Uint size, const Uint offset_l, const Uint offset_r )
    {
      ValueType result = ValueType();
      const ValueType* const lptr = lexpr.data() + offset_l;
      for ( Uint i = 0, j = offset_r; i < size; ++i, ++j )
        result += lptr[i]*rexpr[j];
      return result;
    }
};

// Partial specialization in case the right operand has direct access to the data
// (i.e. the left operand is an expression, and the right operand is not)

template<typename ValueType,typename LExprRef, typename RExprRef>
class vector_vector_product_impl<ValueType,LExprRef,RExprRef,1,0>
{
  public:
    static inline ValueType apply ( const LExprRef lexpr, const RExprRef rexpr,
                                    const Uint size, const Uint offset_l, const Uint offset_r )
    {
      ValueType result = ValueType();
      const ValueType* const rptr = rexpr.data() + offset_r;
      for ( Uint i = offset_l, j = 0; j < size; ++i, ++j )
        result += lexpr[i]*rptr[j];
      return result;
    }
};

// Partial specialization in case both operands are not expressions (i.e. they can directly
// access memory where they store their data

template<typename ValueType,typename LExprRef, typename RExprRef>
class vector_vector_product_impl<ValueType,LExprRef,RExprRef,0,0>
{
  public:
    static inline ValueType apply ( const LExprRef lexpr, const RExprRef rexpr,
                                    const Uint size, const Uint offset_l, const Uint offset_r )
    {
      ValueType result = ValueType();
      const ValueType* const lptr = lexpr.data() + offset_l;
      const ValueType* const rptr = rexpr.data() + offset_r;
      for ( Uint i = 0; i < size; ++i )
        result += lptr[i]*rptr[i];
      return result;
    }
};

#endif

} // Namespace math

} // Namespace pdekit

#endif
