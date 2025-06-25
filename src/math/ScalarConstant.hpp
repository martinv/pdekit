#ifndef PDEKIT_Math_Scalar_Constant_hpp
#define PDEKIT_Math_Scalar_Constant_hpp

#include "common/PDEKit.hpp"

// class for objects that represent scalars

namespace pdekit
{

namespace math
{

template <typename T>
class Scalar
{
};

template <typename T>
class ScalarConstant : public Scalar<ScalarConstant<T>>
{

  public:
  using tensor_type    = ScalarConstant<T>;
  using composite_type = ScalarConstant<T> const &;
  using value_type     = T;
  enum
  {
    TensorRank = 0
  };
  enum
  {
    Size = 1
  };
  enum
  {
    IsExpression = 0
  };

  // constructor initializes value
  ScalarConstant(T const &v) : s(v)
  {
  }

  // for index operations the scalar is the value of each element
  T operator[](const Uint) const
  {
    return s;
  }

  // scalars have zero as size
  Uint size() const
  {
    return 0;
  }

  private:
  T const &s; // value of the scalar
};

} // Namespace math

} // Namespace pdekit

#endif
